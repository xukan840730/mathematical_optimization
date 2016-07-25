#ifndef _BIT_ARRAY_H_
#define _BIT_ARRAY_H_

static const size_t kBitArrayBitsPerBlock = 64;
static const size_t kBitArrayBlockMask = kBitArrayBitsPerBlock - 1;
static const size_t kBitArrayBlockShift = 6;	// val >> 6 == div by 64

#define BIT_ARRAY_ASSERT(a) ASSERT(a)

inline U32 CountSetBitsU32(U32 v)
{
	// bit counting method from http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
	// (in public domain)
	U32 c = v - ((v >> 1) & 0x55555555U);
	c = ((c >> 2) & 0x33333333U) + (c & 0x33333333U);
	c = ((c + (c >> 4) & 0xF0F0F0FU) * 0x1010101U) >> 24;
	return c;
}

inline U64 CountSetBitsU64(U64 v)
{
	// bit counting method from http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
	// (in public domain)
	U64 c = v - ((v >> 1) & 0x5555555555555555ULL);
	c = ((c >>  2) & 0x3333333333333333ULL) + (c & 0x3333333333333333ULL);
	c = ((c >>  4) + c) & 0x0f0f0f0f0f0f0f0fULL;
	c = ((c >>  8) + c) & 0x00ff00ff00ff00ffULL;
	c = ((c >> 16) + c) & 0x0000ffff0000ffffULL;
	c = ((c >> 32) + c) & 0x00000000ffffffffULL;
	return c;
}


inline U32 CountTrailingZeroesU32(U32 value)
{
	U32 ctz = 0;
	if (!value)
		return 32UL;
	if (!(value & 0x0000FFFF))
		ctz += 16, value >>= 16;
	if (!(value & 0x000000FF))
		ctz +=  8, value >>=  8;
	if (!(value & 0x0000000F))
		ctz +=  4, value >>=  4;
	if (!(value & 0x00000003))
		ctz +=  2, value >>=  2;
	if (!(value & 0x00000001))
		ctz +=  1;
	return ctz;
}

inline U64 CountTrailingZeroesU64(U64 value)
{
	U64 ctz = 0;
	if (!value)
		return 64ULL;
	if (!(value & 0x00000000FFFFFFFF))
		ctz += 32, value >>= 32;
	if (!(value & 0x000000000000FFFF))
		ctz += 16, value >>= 16;
	if (!(value & 0x00000000000000FF))
		ctz +=  8, value >>=  8;
	if (!(value & 0x000000000000000F))
		ctz +=  4, value >>=  4;
	if (!(value & 0x0000000000000003))
		ctz +=  2, value >>=  2;
	if (!(value & 0x0000000000000001))
		ctz +=  1;
	return ctz;
}

inline U32 CountLeadingZeroesU32(U32 value)
{
	// binary search count leading 0s
	if (!value)
		return 32UL;

	U32 count = 0;
	for (U32 mask = 0x80000000UL; mask; mask >>= 1)
	{
		if (value & mask)
		{
			return count;
		}
		++count;
	}
	return count;
}

inline U64 CountLeadingZeroesU64(U64 value)
{
	// binary search count leading 0s
	if (!value)
		return 64ULL;

	U32 count = 0;
	for (U64 mask = 0x8000000000000000ULL; mask; mask >>= 1)
	{
		if (value & mask)
		{
			return count;
		}
		++count;
	}
	return count;
}

// -------------------------------------------------------------------------------------------------

template< class STORAGE >
class BitArrayT
{
public:
	typedef STORAGE Storage;

	/* BitArrayT() { } */		// no constructor, we want to use this in unions please!

	void Relocate(ptrdiff_t delta, uintptr_t lowerBound, uintptr_t upperBound)
	{
		// Thanks to a nice feature of C++ templates, a function that is never
		// called will never be instantiated by the compiler. So if your
		// Storage class does not require relocation, it needn't define a
		// Relocate() function -- and any attempt to call Relocate() on the
		// corresponding BitArray will fail to compile (as you'd expect).
		m_storage.Relocate(delta, lowerBound, upperBound);
	}

	inline void SetBlock(size_t iBlock, U64 value) { m_storage.m_block[iBlock] = value; }
	inline U64 GetBlock(size_t iBlock) const { return m_storage.m_block[iBlock]; }

	inline void SetAllBits()				{ AssignAllBits(true); }
	inline void ClearAllBits()				{ AssignAllBits(false); }
	inline void AssignAllBits(bool value)	{ m_storage.Fill((value) ? ~0ULL : 0ULL); }

	void SetBit(U64 index)
	{
		BIT_ARRAY_ASSERT(index < m_storage.GetMaxBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT(index);		// break if watch bit is changed
		m_storage.m_block[index >> kBitArrayBlockShift] |= (1ULL << (index & kBitArrayBlockMask));
	}
	void ClearBit(U64 index)
	{
		BIT_ARRAY_ASSERT(index < m_storage.GetMaxBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT(index);		// break if watch bit is changed
		m_storage.m_block[index >> kBitArrayBlockShift] &= ~(1ULL << (index & kBitArrayBlockMask));
	}

	void AssignBit(U64 index, bool set)
	{
		BIT_ARRAY_ASSERT(index < m_storage.GetMaxBitCount());
		if (set)
			SetBit(index);
		else
			ClearBit(index);
	}

	bool IsBitSet(U64 index) const
	{
		BIT_ARRAY_ASSERT(index < m_storage.GetMaxBitCount());
		return m_storage.m_block[index >> kBitArrayBlockShift] & (1ULL << (index & kBitArrayBlockMask));
	}

	bool operator[](U64 index) const
	{
		return IsBitSet(index);
	}

	size_t GetMaxBitCount() const	{ return m_storage.GetMaxBitCount(); }
	size_t GetNumBlocks() const	{ return m_storage.GetNumBlocks(); }

	size_t CountSetBits() const
	{
		U64 count = 0;
		const size_t numBlocks = GetNumBlocks();
		for (size_t iBlock = 0; iBlock < numBlocks-1; ++iBlock)
		{
			const U64 block = m_storage.m_block[iBlock];
			count += ::CountSetBitsU64(block);
		}

		// special case for final block, which may not be using all 64 bits
		{
			const U64 block = m_storage.m_block[numBlocks-1];
			const U64 mask = BuildBitMask(0, (GetMaxBitCount()-1)%64);
			count += ::CountSetBitsU64(block & mask);
		}

		return (size_t)count;
	}

	bool AreAllBitsClear() const
	{
		const size_t numBlocks = GetNumBlocks();
		for (size_t iBlock = 0; iBlock < numBlocks-1; ++iBlock)
		{
			const U64 block = m_storage.m_block[iBlock];
			if (block != 0ULL)
				return false;
		}

		// special case for final block, which may not be using all 64 bits
		{
			const U64 block = m_storage.m_block[numBlocks-1];
			const U64 mask = BuildBitMask(0, (GetMaxBitCount()-1)%64);
			if ((block & mask) != 0ULL)
				return false;
		}

		return true;
	}

	bool AreAllBitsSet() const
	{
		const size_t numBlocks = GetNumBlocks();
		for (size_t iBlock = 0; iBlock < numBlocks-1; ++iBlock)
		{
			const U64 block = m_storage.m_block[iBlock];
			if (block != ~0ULL)
				return false;
		}

		// special case for final block, which may not be using all 64 bits
		{
			const U64 block = m_storage.m_block[numBlocks-1];
			const U64 mask = BuildBitMask(0, (GetMaxBitCount()-1)%64);
			if ((block & mask) != mask)
				return false;
		}

		return true;
	}

	// Count bits set in range [i0..i1] (inclusive).
	size_t CountBitRangeSetBits(U64 i0, U64 i1) const
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetMaxBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT_RANGE(i0, i1);

		const U32 iBlock0 = i0 >> kBitArrayBlockShift;
		const U32 iBlock1 = i1 >> kBitArrayBlockShift;

		BIT_ARRAY_ASSERT(iBlock0 < m_storage.GetNumBlocks());
		BIT_ARRAY_ASSERT(iBlock1 < m_storage.GetNumBlocks());

		size_t count = 0;

		if (iBlock0 == iBlock1)
		{
			// special case if it's just a single block
			const U64 block = m_storage.m_block[iBlock0];
			const U64 mask = BuildBitMask(i0%64, i1%64);
			count += ::CountSetBitsU64(block & mask);
		}
		else
		{
			for (U32 iBlock = iBlock0+1; iBlock < iBlock1; ++iBlock)
			{
				const U64 block = m_storage.m_block[iBlock];
				count += ::CountSetBitsU64(block);
			}

			{
				const U64 block0 = m_storage.m_block[iBlock0];
				const U64 mask0 = BuildBitMask(i0%64, 63);
				count += ::CountSetBitsU64(block0 & mask0);
			}

			{
				const U64 block1 = m_storage.m_block[iBlock1];
				const U64 mask1 = BuildBitMask(0, i1%64);
				count += ::CountSetBitsU64(block1 & mask1);
			}
		}

		return count;
	}

	void SetBitRange(U64 i0, U64 i1)
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetMaxBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT_RANGE(i0, i1);

		const size_t iBlock0 = i0 >> kBitArrayBlockShift;
		const size_t iBlock1 = i1 >> kBitArrayBlockShift;

		BIT_ARRAY_ASSERT(iBlock0 < m_storage.GetNumBlocks());
		BIT_ARRAY_ASSERT(iBlock1 < m_storage.GetNumBlocks());

		if (iBlock0 == iBlock1)
		{
			// special case if it's just a single block
			const U64 mask = BuildBitMask(i0%64, i1%64);
			m_storage.m_block[iBlock0] |= mask;
		}
		else
		{
			for (size_t iBlock = iBlock0+1; iBlock < iBlock1; ++iBlock)
			{
				m_storage.m_block[iBlock] = ~0ULL;
			}

			{
				const U64 mask0 = BuildBitMask(i0%64, 63);
				m_storage.m_block[iBlock0] |= mask0;
			}

			{
				const U64 mask1 = BuildBitMask(0, i1%64);
				m_storage.m_block[iBlock1] |= mask1;
			}
		}
	}

	void ClearBitRangeUnsafe(U64 i0, U64 i1)
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetAllocatedBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetAllocatedBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT_RANGE(i0, i1);

		const size_t iBlock0 = i0 >> kBitArrayBlockShift;
		const size_t iBlock1 = i1 >> kBitArrayBlockShift;

		BIT_ARRAY_ASSERT(iBlock0 < m_storage.GetNumBlocks());
		BIT_ARRAY_ASSERT(iBlock1 < m_storage.GetNumBlocks());

		if (iBlock0 == iBlock1)
		{
			// special case if it's just a single block
			const U64 mask = BuildBitMask(i0%64, i1%64);
			m_storage.m_block[iBlock0] &= ~mask;
		}
		else
		{
			for (size_t iBlock = iBlock0+1; iBlock < iBlock1; ++iBlock)
			{
				m_storage.m_block[iBlock] = 0ULL;
			}

			{
				const U64 mask0 = BuildBitMask(i0%64, 63);
				m_storage.m_block[iBlock0] &= ~mask0;
			}

			{
				const U64 mask1 = BuildBitMask(0, i1%64);
				m_storage.m_block[iBlock1] &= ~mask1;
			}
		}
	}

	void ClearBitRange(U64 i0, U64 i1)
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetMaxBitCount());
		//BIT_ARRAY_CHECK_WATCH_BIT_RANGE(i0, i1);

		ClearBitRangeUnsafe(i0, i1);
	}

	void ClearUnusedBits()
	{
		ClearBitRangeUnsafe(m_storage.GetAllocatedBitCount()-1, m_storage.GetMaxBitCount());
	}

	bool IsBitRangeAllSet(U64 i0, U64 i1) const
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetMaxBitCount());

		const size_t iBlock0 = i0 >> kBitArrayBlockShift;
		const size_t iBlock1 = i1 >> kBitArrayBlockShift;

		BIT_ARRAY_ASSERT(iBlock0 < m_storage.GetNumBlocks());
		BIT_ARRAY_ASSERT(iBlock1 < m_storage.GetNumBlocks());

		if (iBlock0 == iBlock1)
		{
			// special case if it's just a single block
			const U64 mask = BuildBitMask(i0%64, i1%64);
			return ((m_storage.m_block[iBlock0] & mask) == mask);
		}
		else
		{
			{
				const U64 mask0 = BuildBitMask(i0%64, 63);
				if ((m_storage.m_block[iBlock0] & mask0) != mask0)
					return false;
			}

			for (size_t iBlock = iBlock0+1; iBlock < iBlock1; ++iBlock)
			{
				if (m_storage.m_block[iBlock] != ~0ULL)
					return false;
			}

			{
				const U64 mask1 = BuildBitMask(0, i1%64);
				if ((m_storage.m_block[iBlock1] & mask1) != mask1)
					return false;
			}
		}

		return true;
	}

	bool IsBitRangeAllClearUnsafe(U64 i0, U64 i1) const
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetAllocatedBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetAllocatedBitCount());

		const U32 iBlock0 = i0 >> kBitArrayBlockShift;
		const U32 iBlock1 = i1 >> kBitArrayBlockShift;

		BIT_ARRAY_ASSERT(iBlock0 < m_storage.GetNumBlocks());
		BIT_ARRAY_ASSERT(iBlock1 < m_storage.GetNumBlocks());

		if (iBlock0 == iBlock1)
		{
			// special case if it's just a single block
			const U64 mask = BuildBitMask(i0%64, i1%64);
			return ((m_storage.m_block[iBlock0] & mask) == 0ULL);
		}
		else
		{
			{
				const U64 mask0 = BuildBitMask(i0%64, 63);
				if ((m_storage.m_block[iBlock0] & mask0) != 0ULL)
					return false;
			}

			for (size_t iBlock = iBlock0+1; iBlock < iBlock1; ++iBlock)
			{
				if (m_storage.m_block[iBlock] != 0ULL)
					return false;
			}

			{
				const U64 mask1 = BuildBitMask(0, i1%64);
				if ((m_storage.m_block[iBlock1] & mask1) != 0ULL)
					return false;
			}
		}

		return true;
	}

	bool IsBitRangeAllClear(U64 i0, U64 i1) const
	{
		BIT_ARRAY_ASSERT(i0 <= i1);
		BIT_ARRAY_ASSERT(i0 < m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(i1 < m_storage.GetMaxBitCount());

		return IsBitRangeAllClearUnsafe(i0, i1);
	}

	bool AreUnusedBitsClear() const
	{
		return IsBitRangeAllClearUnsafe(m_storage.GetMaxBitCount(), m_storage.GetAllocatedBitCount()-1);
	}

	U64 FindFirstSetBit() const
	{
		const size_t numBlocks = GetNumBlocks();

		for (size_t i = 0; i < numBlocks - 1; ++i)
		{
			// If the word is non-zero we had set bits... if not we can skip the whole word
			const U64 block = m_storage.m_block[i];

			if (block)
			{
				const U64 lowestBitMask = MaskLowestBitSet(block);
				return GetBitIndexFromMaskedBlock(i, lowestBitMask);
			}
		}

		// special case for final block, which may not be using all 64 bits
		{
			const U64 block = m_storage.m_block[numBlocks-1];

			if (block)
			{
				const U64 mask = BuildBitMask(0, (GetMaxBitCount()-1)%64);
				const U64 maskedBlock = block & mask;

				if (maskedBlock)
				{
					const U64 lowestBitMask = MaskLowestBitSet(maskedBlock);
					return GetBitIndexFromMaskedBlock(numBlocks - 1, lowestBitMask);
				}
			}
		}

		return ~0ULL;
	}

	U64 FindFirstClearBit() const
	{
		const size_t numBlocks = GetNumBlocks();

		for (size_t i = 0; i < numBlocks - 1; ++i)
		{
			// If the inverted word is non-zero we had clear bits... if not we can skip the whole word
			const U64 invBlock = ~m_storage.m_block[i];

			if (invBlock)
			{
				const U64 lowestBitMask = MaskLowestBitSet(invBlock);
				return GetBitIndexFromMaskedBlock(i, lowestBitMask);
			}
		}

		// special case for final block, which may not be using all 64 bits
		{
			const U64 invBlock = ~m_storage.m_block[numBlocks-1];

			if (invBlock)
			{
				const U64 mask = BuildBitMask(0, (GetMaxBitCount()-1)%64);
				const U64 maskedBlock = invBlock & mask;

				if (maskedBlock)
				{
					const U64 lowestBitMask = MaskLowestBitSet(maskedBlock);
					return GetBitIndexFromMaskedBlock(numBlocks - 1, lowestBitMask);
				}
			}
		}

		return ~0ULL;
	}
	
	// ----------------------------------------------------------------------------------------------------
	// FUNCTIONS PAST HERE DON'T NECESSARILY WORK WITH NON-64-ALIGNED SIZES
	// ----------------------------------------------------------------------------------------------------

	U64 FindLastSetBit() const
	{
		const size_t numBlocks = GetNumBlocks();
		for (int i = numBlocks - 1; i >= 0; --i)
		{
			// If the word is non-zero we had set bits... if not we can skip the whole word
			const U64 currentBlock = m_storage.m_block[i];
			if (currentBlock)
			{
				const U64 maskedBlock = MaskHighestBitSet(currentBlock);
				return GetBitIndexFromMaskedBlock(i, maskedBlock);
			}
		}

		return ~0ULL;
	}

	U64 FindNextSetBit(U64 currentIndex) const
	{
		U64 currentBlock = currentIndex >> kBitArrayBlockShift;

		// Set the corresponding bit.
		U64 currentBitMask = 1ULL << (currentIndex & kBitArrayBlockMask);

		// Fill in all bits lower than the bit just set.
		U64 currentBitMaskFilled = currentBitMask | (currentBitMask - 1ULL);

		// Negate and create a mask to filter out all higher bits.
		U64 currentBlockMask = ~currentBitMaskFilled;

		// Mask the current word and then test it.
		U64 maskedBlock = MaskLowestBitSet(m_storage.m_block[currentBlock] & currentBlockMask);
		if (maskedBlock)
		{
			return GetBitIndexFromMaskedBlock(currentBlock, maskedBlock);
		}

		currentBlock++;
		size_t numBlocks = GetNumBlocks();
		for (; currentBlock < numBlocks; ++currentBlock)
		{
			// If the word is non-zero we had set bits... if not we can skip the whole word
			maskedBlock = MaskLowestBitSet(m_storage.m_block[currentBlock]);
			if (maskedBlock)
			{
				return GetBitIndexFromMaskedBlock(currentBlock, maskedBlock);
			}
		}

		return ~0ULL;
	}

	U64 FindNextClearBit(U64 currentIndex) const
	{
		U64 currentBlock = currentIndex >> kBitArrayBlockShift;

		// Set the corresponding bit.
		U64 currentBitMask = 1ULL << (currentIndex & kBitArrayBlockMask);

		// Fill in all bits lower than the bit just set.
		U64 currentBitMaskFilled = currentBitMask | (currentBitMask - 1ULL);

		// Negate and create a mask to filter out all higher bits.
		U64 currentBlockMask = ~currentBitMaskFilled;

		// Mask the current word and then test it.
		U64 maskedBlock = MaskLowestBitSet((~m_storage.m_block[currentBlock]) & currentBlockMask);
		if (maskedBlock)
		{
			return GetBitIndexFromMaskedBlock(currentBlock, maskedBlock);
		}

		currentBlock++;
		size_t numBlocks = GetNumBlocks();
		for (; currentBlock < numBlocks; ++currentBlock)
		{
			// If the word is non-zero we had set bits... if not we can skip the whole word
			maskedBlock = MaskLowestBitSet(~m_storage.m_block[currentBlock]);
			if (maskedBlock)
			{
				return GetBitIndexFromMaskedBlock(currentBlock, maskedBlock);
			}
		}

		return ~0ULL;
	}
	
	//
	// Find the first bit such that it and the next `numBits` bits are clear.
	//
	//U64 FindFirstNClearBitsByBlock(U32 numBits, U32 startBlock, U32 endBlock) const
	//{
	//	U32 foundBits = 0;
	//		
	//	for (U32 iBlock = startBlock; iBlock <= endBlock; iBlock++)
	//	{
	//		U64 x = m_storage.m_block[iBlock];

	//		if (x == 0ULL)
	//		{
	//			if (foundBits + kBitArrayBitsPerBlock >= numBits)
	//				return iBlock * kBitArrayBitsPerBlock - foundBits;
	//			else
	//				foundBits += kBitArrayBitsPerBlock;
	//		}
	//		else
	//		{
	//			// We have a bit set in here.
	//			{
	//				U32 trailingBits = CountTrailingZeroesU64(x);
//
//					if (foundBits + trailingBits >= numBits)
//						return iBlock * kBitArrayBitsPerBlock - foundBits;
//					else
//						foundBits = 0;
//				}
//
//				// Otherwise, if numBits < 64, find n zeros.
	//			if (numBits < kBitArrayBitsPerBlock)
	//			{
	//				U32 trailingBits = FindTrailingNClearBitsU64(x, numBits);
	//				if (trailingBits < kBitArrayBitsPerBlock)
	//				{
	//					return (iBlock * kBitArrayBitsPerBlock + trailingBits);
	//				}
	//			}
//
	//			foundBits = CountLeadingZeroesU64(x);
	//			if (foundBits >= numBits)
	//				return (iBlock * kBitArrayBitsPerBlock + kBitArrayBitsPerBlock - foundBits);
	//		}
	//	}
//
//		return ~0ULL;
	//}
	
//	U64 FindFirstNClearBits(U32 numBits) const
//	{
//		return FindFirstNClearBitsByBlock(numBits, 0, GetNumBlocks() - 1);
//	}

	//
	// Find the first bit such that it and the next `numBits` bits are clear.
	//
//	U64 FindFirstNClearBitsByBlockAligned(U32 numBits, U32 startBlock, U32 endBlock, U64 alignment) const
//	{
//		U32 foundBits = 0;
//
//		U64 blockAlignment = alignment / 64; // 0 => All blocks. 1 => iBlock & 1 == 0, 2 => iBlock & 3 == 0, etc.
//
//		U64 validBits =
//			  alignment == 32 ? 0x0000000100000001ull
//			: alignment == 16 ? 0x0001000100010001ull
//			: alignment ==  8 ? 0x0101010101010101ull
//			: alignment ==  4 ? 0x1111111111111111ull
//			: alignment ==  2 ? 0x5555555555555555ull
//			:                   0xFFFFFFFFFFFFFFFFull;
//
//		
//		for (U32 iBlock = startBlock; iBlock <= endBlock; iBlock++)
//		{
//			if (!foundBits && blockAlignment && ((iBlock & ((1 << blockAlignment)-1)) != 0))
//				continue; // Can't start on this block
//
//			U64 x = m_storage.m_block[iBlock];
//			
//			if (x == 0ULL)
//			{
//				if (foundBits + kBitArrayBitsPerBlock >= numBits)
//					return iBlock * kBitArrayBitsPerBlock - foundBits;
//				else
//					foundBits += kBitArrayBitsPerBlock;
//			}
//			else
//			{
//				// We have a bit set in here.
//
//				if (foundBits) // Restricted since otherwise alignment won't take effect
//				{
//					// Bits in the low part (bit 0, 1, ... )					
//					U32 trailingBits = CountTrailingZeroesU64(x);
//
//					if (foundBits + trailingBits >= numBits)
//						return iBlock * kBitArrayBitsPerBlock - foundBits;
//					else
//						foundBits = 0;
//				}
//
//				// Otherwise, if numBits < 64, find n zeros.
//				if (numBits < kBitArrayBitsPerBlock)
//				{
//					U32 trailingBits = FindTrailingNClearBitsU64(x, numBits, validBits);
//					if (trailingBits < kBitArrayBitsPerBlock)
//					{
//						return (iBlock * kBitArrayBitsPerBlock + trailingBits);
//					}
//				}
//
//				// Bits in the high part (bit 63, 62, ... )
//				U32 leadingBits = CountLeadingZeroesU64(x);
//				// Only allow number according to alignment
//				foundBits = ALIGN_SIZE_LOWER(leadingBits, alignment);
//				// This should have been cought above, but I guess it doesn't hurt
//				if (foundBits >= numBits)
//					return (iBlock * kBitArrayBitsPerBlock + kBitArrayBitsPerBlock - foundBits);
//			}
//		}
//
//		return ~0ULL;
//	}
	
//	U64 FindFirstNClearBitsAligned(U32 numBits, U32 alignment) const
//	{
//		return FindFirstNClearBitsByBlockAligned(numBits, 0, GetNumBlocks() - 1, alignment);
//	}	

	//
	// Find the first bit such that it and the next `numBits` bits are set.
	//	
//	U64 FindFirstNSetBitsByBlock(U32 numBits, U32 startBlock, U32 endBlock) const
//	{
//		U32 foundBits = 0;
//			
//		for (U32 iBlock = startBlock; iBlock <= endBlock; iBlock++)
//		{
//			U64 x = m_storage.m_block[iBlock];
//
//			if (x == ~0ULL)
//			{
//				if (foundBits + kBitArrayBitsPerBlock >= numBits)
//					return iBlock * kBitArrayBitsPerBlock - foundBits;
//				else
//					foundBits += kBitArrayBitsPerBlock;
//			}
//			else
//			{
//				// We have a bit cleared in here.
//				// If foundBits > 0, then we must find enough trailing ones.
//				if (foundBits > 0)
//				{
//					U32 trailingBits = CountTrailingZeroesU64(~x);
///
//					if (foundBits + trailingBits >= numBits)
//						return iBlock * kBitArrayBitsPerBlock - foundBits;
//					else
//						foundBits = 0;
//				}
//
//				// Otherwise, if numBits < 64, find n ones.
//				if (numBits < kBitArrayBitsPerBlock)
//				{
//					U32 trailingBits = FindTrailingNSetBitsU64(x, numBits);
//					if (trailingBits < kBitArrayBitsPerBlock)
//					{
//						return (iBlock * kBitArrayBitsPerBlock + trailingBits);
//					}
//				}
//
//				foundBits = CountLeadingZeroesU64(~x);
//				if (foundBits >= numBits)
//					return (iBlock * kBitArrayBitsPerBlock + kBitArrayBitsPerBlock - foundBits);
//			}
//		}
//
//		return ~0ULL;
//	}

//	U64 FindFirstNSetBits(U32 numBits) const
//	{
//		return FindFirstNSetBitsByBlock(numBits, 0, GetNumBlocks() - 1);
//	}
//	
//	size_t CountLeadingZeroes() const
//	{
//		size_t numLeadingZeroes = 0;
//		const size_t numBlocks = GetNumBlocks();
//		for (int i = numBlocks - 1; i >= 0; --i)
//		{
//			const U64 currentBlock = m_storage.m_block[i];
//			if (currentBlock == 0ULL)
//			{
//				numLeadingZeroes += kBitArrayBitsPerBlock; // all 64 bits were zeroes
//			}
//			else
//			{
//				numLeadingZeroes += CountLeadingZeroesU64(currentBlock);
//				break;
//			}
//		}
//
//		return numLeadingZeroes;
//	}

	static void BitwiseOr(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		const size_t numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		const size_t repeatCount = numBlocks / 4;
		const size_t finalOffset = repeatCount * 4;

		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);

		// unrolling this to see if we can improve SPU performance
		for (size_t i = 0; i < repeatCount * 4; i += 4)
		{
			const U64 block0 = src1.m_storage.m_block[i+0] | src2.m_storage.m_block[i+0];
			const U64 block1 = src1.m_storage.m_block[i+1] | src2.m_storage.m_block[i+1];
			const U64 block2 = src1.m_storage.m_block[i+2] | src2.m_storage.m_block[i+2];
			const U64 block3 = src1.m_storage.m_block[i+3] | src2.m_storage.m_block[i+3];

			pDest->m_storage.m_block[i+0] = block0;
			pDest->m_storage.m_block[i+1] = block1;
			pDest->m_storage.m_block[i+2] = block2;
			pDest->m_storage.m_block[i+3] = block3;
		}

		switch (numBlocks % 4)
		{
			case 3:	pDest->m_storage.m_block[finalOffset + 2] = src1.m_storage.m_block[finalOffset + 2] | src2.m_storage.m_block[finalOffset + 2];
			case 2:	pDest->m_storage.m_block[finalOffset + 1] = src1.m_storage.m_block[finalOffset + 1] | src2.m_storage.m_block[finalOffset + 1];
			case 1:	pDest->m_storage.m_block[finalOffset + 0] = src1.m_storage.m_block[finalOffset + 0] | src2.m_storage.m_block[finalOffset + 0];
			case 0:	break;
		}
	}

	static void BitwiseXor(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		size_t numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = src1.m_storage.m_block[i] ^ src2.m_storage.m_block[i];
		}
	}

	static void BitwiseAnd(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		size_t numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = src1.m_storage.m_block[i] & src2.m_storage.m_block[i];
		}
	}

	static void BitwiseAndComp(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		size_t numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = src1.m_storage.m_block[i] & ~src2.m_storage.m_block[i];
		}
	}

	static void BitwiseOrComp(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		size_t numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = src1.m_storage.m_block[i] | ~src2.m_storage.m_block[i];
		}
	}

	static void BitwiseCompAndComp(BitArrayT* pDest, const BitArrayT& src1, const BitArrayT& src2)
	{
		U32 numBlocks = Min(src1.GetNumBlocks(), src2.GetNumBlocks());
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for(U32 i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = ~src1.m_storage.m_block[i] & ~src2.m_storage.m_block[i];
		}
	}

	static void BitwiseNot(BitArrayT* pDest, const BitArrayT& src)
	{
		size_t numBlocks = src.GetNumBlocks();
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = ~src.m_storage.m_block[i];
		}
	}

	static void Copy(BitArrayT* pDest, const BitArrayT& other)
	{
		size_t numBlocks = pDest->GetNumBlocks();
		BIT_ARRAY_ASSERT(numBlocks >= other.GetNumBlocks());
		size_t minNumBlocks = Min(numBlocks, other.GetNumBlocks());
		
		size_t i = 0;
		for (; i < minNumBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = other.m_storage.m_block[i];
		}

		for (; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] = 0;
		}
	}

	static void BitwiseSelect(BitArrayT* pDest, const BitArrayT& srcIf0, const BitArrayT& srcIf1, const BitArrayT& srcSelect)
	{
		const size_t numBlocks = srcSelect.GetNumBlocks();
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);
		BIT_ARRAY_ASSERT(srcIf0.GetNumBlocks() >= numBlocks);
		BIT_ARRAY_ASSERT(srcIf1.GetNumBlocks() >= numBlocks);
		for (size_t i = 0; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i] =
				(srcIf0.m_storage.m_block[i] & (~srcSelect.m_storage.m_block[i])) |
				(srcIf1.m_storage.m_block[i] &  (srcSelect.m_storage.m_block[i]));
		}
	}

	static void AddWithCarry(BitArrayT* pDest, const BitArrayT& src, U64 addend)
	{
		size_t numBlocks = src.GetNumBlocks();
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);

		U64 carry = addend;
		for (size_t i = 0; i < numBlocks; ++i)
		{
			const U64 result = src.m_storage.m_block[i] + carry;
			carry = (result < carry); // carry if necessary
			pDest->m_storage.m_block[i] = result;
		}
	}

	static void SubtractWithBorrow(BitArrayT* pDest, const BitArrayT& src, U64 subtrahend)
	{
		size_t numBlocks = src.GetNumBlocks();
		BIT_ARRAY_ASSERT(pDest);
		BIT_ARRAY_ASSERT(pDest->GetNumBlocks() >= numBlocks);

		U64 borrow = subtrahend;
		for (size_t i = 0; i < numBlocks; ++i)
		{
			const U64 block = src.m_storage.m_block[i];
			const U64 result = block - borrow; // no need to actually borrow -- this works because we're using unsigned arithmetic
			borrow = (block < borrow); // borrow if necessary
			pDest->m_storage.m_block[i] = result;
		}
	}

	static bool IsEqual(const BitArrayT& src1, const BitArrayT& src2)
	{
		U64 temp = 0;
		size_t numBlocks1 = src1.GetNumBlocks();
		size_t numBlocks2 = src2.GetNumBlocks();

		if (numBlocks1 != numBlocks2)
		{
			return false;
		}

		for (size_t iBlock = 0; iBlock < numBlocks1; ++iBlock)
		{
			const U64 vMask = (iBlock == numBlocks1 - 1) ? BuildBitMask(0, (src1.GetMaxBitCount() - 1) % 64) : ~0ULL;
			temp |= vMask & (src1.m_storage.m_block[iBlock] ^ src2.m_storage.m_block[iBlock]);
		}

		return temp == 0;
	}

	static void ShortShiftLeft(BitArrayT* pDest, const BitArrayT& src, size_t bitShift)
	{
		const size_t numBlocks = src.GetNumBlocks();

		BIT_ARRAY_ASSERT(pDest != &src);		// can't do this in place, sorry kids
		BIT_ARRAY_ASSERT(bitShift < kBitArrayBitsPerBlock);
		BIT_ARRAY_ASSERT(bitShift < src.m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(pDest->m_storage.GetMaxBitCount() == src.m_storage.GetMaxBitCount());

		pDest->m_storage.m_block[0] = src.m_storage.m_block[0] << bitShift;

		for (size_t i = 1; i < numBlocks; ++i)
		{
			pDest->m_storage.m_block[i]  = src.m_storage.m_block[i]   << bitShift;
			pDest->m_storage.m_block[i] |= src.m_storage.m_block[i-1] >> (kBitArrayBitsPerBlock - bitShift);
		}

		// the most significant bits of the result will be garbage, so let's
		// clear them in case we do something with them
		pDest->ClearUnusedBits();
	}

	static void ShortShiftRight(BitArrayT* pDest, const BitArrayT& src, size_t bitShift)
	{
		const U32 numBlocks = src.GetNumBlocks();

		BIT_ARRAY_ASSERT(pDest != &src);		// can't do this in place, sorry kids
		BIT_ARRAY_ASSERT(bitShift < kBitArrayBitsPerBlock);
		BIT_ARRAY_ASSERT(bitShift < src.m_storage.GetMaxBitCount());
		BIT_ARRAY_ASSERT(pDest->m_storage.GetMaxBitCount() == src.m_storage.GetMaxBitCount());

		for (U32 i = 0; i < numBlocks-1; ++i)
		{
			pDest->m_storage.m_block[i]  = src.m_storage.m_block[i] >> bitShift;
			pDest->m_storage.m_block[i] |= src.m_storage.m_block[i+1] << (kBitArrayBitsPerBlock - bitShift);
		}

		pDest->m_storage.m_block[numBlocks-1] = src.m_storage.m_block[numBlocks-1] >> bitShift;
	}

	// NB: This is NOT a "robust" iterator; you should not modify the bit array
	// while iterating through it. (Because we cache blocks, there's a chance
	// we'll miss the bit you modify).
	class Iterator
	{
		typedef BitArrayT<STORAGE> Bits;

		const Bits*	m_pBits;
		U64			m_currentBit;
		U32			m_numBlocks;
		U32			m_currentBlockIndex;
		U64			m_currentBlock;			// this has the bits removed from it one by one

	public:
		Iterator(const Bits& bits) : m_pBits(&bits)
		{
			m_numBlocks = m_pBits->GetNumBlocks();
			m_currentBit = End();
		}

		Iterator(const Bits* pBits) : m_pBits(pBits)
		{
			m_numBlocks = m_pBits->GetNumBlocks();
			m_currentBit = End();
		}

		Iterator()
		{
			m_currentBit = ~0;
			m_numBlocks = 0;
			m_pBits = NULL;
		}

		void Relocate(ptrdiff_t deltaPos, uintptr_t lowerBound, uintptr_t upperBound)
		{
			RelocatePointer(m_pBits, deltaPos, lowerBound, upperBound);
		}

		Iterator EndIter() const
		{
			Iterator endIter = *this;
			endIter.m_currentBit = End();
			return endIter;
		}

		U64 First()
		{
			BIT_ARRAY_ASSERT(m_pBits);
			for (m_currentBlockIndex = 0; m_currentBlockIndex < m_numBlocks; ++m_currentBlockIndex)
			{
				m_currentBlock = m_pBits->GetBlock(m_currentBlockIndex);

				if (m_currentBlock)
				{
					const U64 currentBitMask = m_pBits->MaskLowestBitSet(m_currentBlock);
					m_currentBlock ^= currentBitMask;
					m_currentBit = m_pBits->GetBitIndexFromMaskedBlock(m_currentBlockIndex, currentBitMask);
					break;
				}
			}

			return m_currentBit;
		}

		U64 Current() const
		{
			BIT_ARRAY_ASSERT(m_pBits);
			return m_currentBit;
		}

		U64 End() const
		{
			BIT_ARRAY_ASSERT(m_pBits);
			return m_pBits->GetMaxBitCount();
		}

		U64 Advance()
		{
			BIT_ARRAY_ASSERT(m_pBits);
			if (m_currentBlock)
			{
				const U64 currentBitMask = m_pBits->MaskLowestBitSet(m_currentBlock);
				m_currentBlock ^= currentBitMask;
				m_currentBit = m_pBits->GetBitIndexFromMaskedBlock(m_currentBlockIndex, currentBitMask);
				return m_currentBit;
			}

			m_currentBlockIndex++;

			for (; m_currentBlockIndex < m_numBlocks; ++m_currentBlockIndex)
			{
				m_currentBlock = m_pBits->GetBlock(m_currentBlockIndex);

				if (m_currentBlock)
				{
					const U64 currentBitMask = m_pBits->MaskLowestBitSet(m_currentBlock);
					m_currentBlock ^= currentBitMask;
					m_currentBit = m_pBits->GetBitIndexFromMaskedBlock(m_currentBlockIndex, currentBitMask);
					return m_currentBit;
				}
			}

			m_currentBit = End();
			return m_currentBit;
		}
	};

	static U64 BuildBitMask(U64 i0, U64 i1)
	{
		// this relies on left/right shift operators filling with 0
		return (~0ULL<<i0) & (~0ULL>>(63-i1));
	}
	const Storage& GetStorage() const { return m_storage; }
protected:

	U64 MaskLowestBitSet(U64 word) const
	{
#if defined NDI_PLAT_ORBIS
		return __blsi_u64(word);
#else
		return word & (-word);
#endif
	}

	U64 MaskHighestBitSet(U64 word) const
	{
		if (word != 0ULL)
			return (1ULL << ((kBitArrayBitsPerBlock - 1) - (int)CountLeadingZeroesU64(word)));
		else
			return 0ULL;
	}

	U64 GetBitIndexFromMaskedBlock(U64 blockIndex, U64 maskedBlock) const
	{
		return blockIndex * kBitArrayBitsPerBlock + kBitArrayBitsPerBlock - 1ULL - CountLeadingZeroesU64(maskedBlock);
	}

	Storage m_storage;
};

// -------------------------------------------------------------------------------------------------

// Implements a BitArrayT type that stores its data externally.
struct ExternalBitArrayStorage
{
	static const size_t DetermineNumBlocks(size_t maxBits)
	{
		// support any input value -- round up to nearest multiple of block size
		return (maxBits + kBitArrayBlockMask) / kBitArrayBitsPerBlock;
	}

	static const size_t DetermineCapacity(size_t maxBits)
	{
		return DetermineNumBlocks(maxBits) * kBitArrayBitsPerBlock;
	}

	ExternalBitArrayStorage()
	{
		m_block = NULL;
		m_maxBits = m_numBlocks = 0U;
		m_allocatedBits = 0;

		#if BIT_ARRAY_DEBUG
			m_watchIndex = ~0U;
		#endif
	}

	void Init(U64 maxBits, U64* aBlock)
	{
		BIT_ARRAY_ASSERT(maxBits > 0);
		BIT_ARRAY_ASSERT(aBlock != NULL);

		m_maxBits = maxBits;
		m_numBlocks = DetermineNumBlocks(maxBits);
		m_block = aBlock;
		m_allocatedBits = DetermineCapacity(maxBits);

		BIT_ARRAY_ASSERT(((U64)m_allocatedBits & kBitArrayBlockMask) == 0);
	}

	//void Relocate(ptrdiff_t delta, uintptr_t lowerBound, uintptr_t upperBound)
	//{
	//	RelocatePointer(m_block, delta, lowerBound, upperBound);
	//}

	size_t GetMaxBitCount() const { return m_maxBits; }
	size_t GetAllocatedBitCount() const { return m_allocatedBits; }

	size_t GetNumBlocks() const { return m_numBlocks; }

	void Fill(U64 value) { memset(m_block, value, sizeof(U64)*m_numBlocks); } 

	U64*	m_block;			// U64 m_block[m_numBlocks];
	U32		m_numBlocks;
	U64		m_maxBits;
	U64		m_allocatedBits;

#if BIT_ARRAY_DEBUG
	U32		m_watchIndex;
#endif
};


// A BitArrayT type that stores its data externally.
class ExternalBitArray : public BitArrayT< ExternalBitArrayStorage >
{
public:
	typedef BitArrayT< ExternalBitArrayStorage > ParentClass;

	// an array of this many U64's should be passed to Init()
	static const size_t DetermineNumBlocks(size_t maxBits)
	{
		return Storage::DetermineNumBlocks(maxBits);
	}

	// round up to nearest multiple of block size
	static const size_t DetermineCapacity(size_t maxBits)
	{
		return Storage::DetermineCapacity(maxBits);
	}

	ExternalBitArray() : ParentClass()
	{
		//AssignAllBits(false); // no need -- we have zero bits until Init() is called
	}

	ExternalBitArray(U64 maxBits, U64* aBlock, bool flag = false)
	{
		m_storage.Init(maxBits, aBlock);
		AssignAllBits(flag);
	}

	void Init(U64 maxBits, U64* aBlock, bool flag = false)
	{
		m_storage.Init(maxBits, aBlock);
		AssignAllBits(flag);
	}

	// used when copy bit-array, can't set/clear all bits.
	void InitNoAssign(U64 maxBits, U64* aBlock)
	{
		m_storage.Init(maxBits, aBlock);
	}
};


#endif
