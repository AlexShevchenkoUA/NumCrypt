#include"Interface.h"

//class BigInteger

void copy(const BigInteger& scr, BigInteger& dst) {
	for (int i = 0; i < dst.capacity; i++)
		dst.number[i] = scr.number[i];
}

void copy(const BigDouble& scr, BigInteger& dst) {
	for (int i = 0; i < dst.capacity; i++)
		dst.number[i] = scr.number[i];
}

void copy(const BigInteger& scr, BigDouble& dst) {
	dst.nullify();
	for (int i = 0; i < scr.capacity; i++)
		dst.number[i] = scr.number[i];
}

int highestNonZeroBit(const BigInteger& x) {
	int i, j;
	for (i = x.capacity - 1; i > -1; i--)
		if (x.number[i] != 0)
			break;
	for (j = x.base - 1; j > -1; j--)
		if (((x.number[i] >> j) & 1) != 0)
			break;
	return (x.base * i + j);
}

int highestNonZeroBlock(const BigInteger& x) {
	int i;
	for (i = x.capacity - 1; i > 0; i--)
		if (x.number[i] != 0)
			break;
	return i;
}

void add(const BigInteger& x, const BigInteger& y, BigInteger& r, unsigned int& carry) {
	r.nullify();
	carry = 0;
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] + y.number[i] + carry;
		carry = ((x.number[i] & y.number[i]) | ((x.number[i] | y.number[i]) & (~r.number[i]))) >> 31;
	}
}

void add(const BigInteger& x, const BigInteger& y, BigInteger& r) {
	r.nullify();
	unsigned int carry = 0;
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] + y.number[i] + carry;
		carry = ((x.number[i] & y.number[i]) | ((x.number[i] | y.number[i]) & (~r.number[i]))) >> 31;
	}
}

void add(const BigInteger& x, const BigInteger& y, BigDouble& r) {
	r.nullify();
	unsigned int carry = 0;
	for (int i = 0; i < x.capacity; i++) {
		r.number[i] = x.number[i] + y.number[i] + carry;
		carry = ((x.number[i] & y.number[i]) | ((x.number[i] | y.number[i]) & (~r.number[i]))) >> 31;
	}
	r.number[x.capacity] = carry;
}

void sub(const BigInteger& x, const BigInteger& y, BigInteger& r, unsigned int& borrow) {
	borrow = 0;
	r.nullify();
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] - y.number[i] - borrow;
		borrow = (((~x.number[i]) & y.number[i]) | (((~x.number[i]) | y.number[i]) & r.number[i])) >> 31;
	}
}

void sub(const BigInteger& x, const BigInteger& y, BigInteger& r) {
	r.nullify();
	unsigned int borrow = 0;
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] - y.number[i] - borrow;
		borrow = (((~x.number[i]) & y.number[i]) | (((~x.number[i]) | y.number[i]) & r.number[i])) >> 31;
	}
}

void increment(const BigInteger& x, BigInteger& r) {
	copy(x, r);
	unsigned int carry = 1;
	int i = 0;
	while ((carry == 1) && (i < x.capacity)) {
		r.number[i]++;
		if (r.number[i] != 0)
			carry = 0;
		i++;
	}
}

void mult(const BigInteger& x, const unsigned int y, BigInteger& r) {
	r.nullify();
	unsigned int hh, hl, lh, ll;
	unsigned int temp;
	unsigned int overflow = 0;
	for (int i = 0; i < r.capacity; i++) {
		hh = (x.number[i] >> 16) * (y >> 16);
		hl = (x.number[i] >> 16) * (y & 0xFFFF);
		lh = (x.number[i] & 0xFFFF) * (y >> 16);
		ll = (x.number[i] & 0xFFFF) * (y & 0xFFFF);
		r.number[i] = ll + ((hl + lh) << 16) + overflow;
		temp = overflow;
		overflow = hh + (hl >> 16) + (lh >> 16);
		//for another part of cycle use unsigned int variable "hh" as a buffer for variable "overflow"
		hh = temp;
		temp = (hl & 0xFFFF) + (lh & 0xFFFF);
		overflow = overflow + (temp >> 16);
		temp = (hl + lh) << 16;
		overflow = overflow + (((ll & temp) | ((ll | temp) & (~(ll + temp)))) >> 31);
		temp = temp + ll;
		overflow = overflow + (((hh & temp) | ((hh | temp) & (~(hh + temp)))) >> 31);
	}
}

void mult(const BigInteger& x, const BigInteger& y, BigInteger& r) {
	r.nullify();
	BigInteger temp(0);
	BigInteger aux(0);
	for (int i = 0; i < y.capacity; i++) {
		mult(x, y.number[i], temp);
		leftBlocksShift(temp, i, aux);
		copy(r, temp);
		add(temp, aux, r);
	}
}

void mult(const BigInteger& x, const BigInteger& y, BigDouble& r) {
	r.nullify();
	BigDouble temp_x;
	copy(x, temp_x);
	BigDouble temp(0);
	BigDouble aux(0);
	for (int i = 0; i < y.capacity; i++) {
		mult(temp_x, y.number[i], temp);
		leftBlocksShift(temp, i, aux);
		copy(r, temp);
		add(temp, aux, r);
	}
}

void div(const BigInteger& x, const BigInteger& y, BigInteger& q, BigInteger& r) {
	BigInteger temp(0);
	BigInteger aux(0);
	int k = highestNonZeroBit(y) + 1;
	copy(x, r);
	q.nullify();
	while (cmp(r, y) != -1) {
		int t = r.highestNonZeroBit() + 1;
		leftShift(y, t - k, temp);
		if (cmp(r, temp) == -1) {
			t--;
			leftShift(y, t - k, temp);
		}
		copy(r, aux);
		sub(aux, temp, r);
		q.setBit(t - k);
	}
}

void mod(const BigInteger& x, const BigInteger& n, BigInteger& r) {
	BigInteger temp(0);
	BigInteger aux(0);
	int k = highestNonZeroBit(n) + 1;
	copy(x, r);
	while (cmp(r, n) != -1) {
		int t = r.highestNonZeroBit() + 1;
		leftShift(n, t - k, temp);
		if (cmp(r, temp) == -1) {
			t--;
			leftShift(n, t - k, temp);
		}
		copy(r, aux);
		sub(aux, temp, r);
	}
}

void power(const BigInteger& x, const BigInteger& e, BigInteger& r) {
	BigInteger temp(0);
	r.set(1);
	int blocks = highestNonZeroBlock(e);
	for (int i = blocks; i > -1; i--) {
		for (int j = e.base - 1; j > -1; j--) {
			if (((e.number[i] >> j) & 0x1) == 1) {
				copy(r, temp);
				mult(x, temp, r);
			}
			if ((i + j) != 0) {
				copy(r, temp);
				mult(temp, temp, r);
			}
		}
	}
}

void powerWindow(const BigInteger& x, const BigInteger& e, BigInteger& r) {
	BigInteger temp(1);
	r.set(1);
	BigInteger powers[16];
	powers[0].set(1);
	copy(x, powers[1]);
	for (int i = 2; i < 16; i++)
		mult(x, powers[i - 1], powers[i]);
	int blocks = highestNonZeroBlock(e);
	for (int i = blocks; i > -1; i--) {
		for (int j = (e.base / 4) - 1; j > -1; j--) {
			copy(r, temp);
			mult(powers[(e.number[i] >> (4 * j)) & 0xF], temp, r);
			if ((i + j) != 0) {
				for (int k = 0; k < 4; k++) {
					copy(r, temp);
					mult(temp, temp, r);
				}
			}
		}
	}
}

int cmp(const BigInteger& x, const BigInteger& y) {
	int i = x.capacity - 1;
	int res;
	while ((x.number[i] == y.number[i]) && i > -1)
		i--;
	if (i == -1)
		res = 0;
	else if (x.number[i] > y.number[i])
		res = 1;
	else
		res = -1;
	return res;
}

void leftShift(const BigInteger& x, int n, BigInteger& r) {
	r.nullify();
	if (n >= (x.capacity * x.base))
		return;
	int sub_block = n % x.base;
	int blocks = (n - sub_block) / x.base;
	int left_lim = x.capacity - blocks;
	unsigned int param;
	if (sub_block == 0)
		param = 0;
	else
		param = 0xFFFFFFFF;
	for (int i = 1; i < left_lim; i++) {
		r.number[x.capacity - i] = (x.number[left_lim - i] << sub_block) +
			(param & (x.number[left_lim - (i + 1)] >> (x.base - sub_block)));
	}
	r.number[blocks] = x.number[0] << sub_block;
}

void leftBlocksShift(const BigInteger& x, int n, BigInteger& r) {
	r.nullify();
	if (n >= x.capacity)
		return;
	for (int i = x.capacity - 1; i >= n; i--)
		r.number[i] = x.number[i - n];
}

void rightShift(const BigInteger& x, int n, BigInteger& r) {
	r.nullify();
	if (n >= (x.base * x.capacity))
		return;
	int sub_block = n % x.base;
	int blocks = n / x.base;
	int right_lim = x.capacity - blocks;
	unsigned int param;
	if (sub_block == 0)
		param = 0;
	else
		param = 0xFFFFFFFF >> (x.base - sub_block);
	for (int i = 0; i < right_lim - 1; i++)
		r.number[i] = (x.number[blocks + i] >> sub_block) + ((param & (x.number[blocks + i + 1])) << (x.base - sub_block));
	r.number[right_lim - 1] = x.number[x.capacity - 1] >> sub_block;
}

void rightBlocksShift(const BigInteger& x, int n, BigInteger& r) {
	r.nullify();
	if (n >= x.capacity) 
		return;
	for (int i = 0; i < x.capacity - n; i++)
		r.number[i] = x.number[i + n];
}

//class BigDouble

void copy(const BigDouble& scr, BigDouble& dst) {
	for (int i = 0; i < dst.capacity; i++)
		dst.number[i] = scr.number[i];
}

int highestNonZeroBit(const BigDouble& x) {
	int i, j;
	for (i = x.capacity - 1; i > -1; i--) {
		if (x.number[i] != 0)
			break;
	}
	for (j = x.base - 1; j > -1; j--) {
		if (((x.number[i] >> j) & 1) != 0)
			break;
	}
	return (x.base * i + j);
}

int highestNonZeroBlock(const BigDouble& x) {
	int i;
	for (i = x.capacity - 1; i > 0; i--)
		if (x.number[i] != 0)
			break;
	return i;
}

void add(const BigDouble& x, const BigDouble& y, BigDouble& r) {
	r.nullify();
	unsigned int carry = 0;
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] + y.number[i] + carry;
		carry = ((x.number[i] & y.number[i]) | ((x.number[i] | y.number[i]) & (~r.number[i]))) >> 31;
	}
}

void sub(const BigDouble& x, const BigDouble& y, BigDouble& r) {
	unsigned int borrow = 0;
	r.nullify();
	for (int i = 0; i < r.capacity; i++) {
		r.number[i] = x.number[i] - y.number[i] - borrow;
		borrow = (((~x.number[i]) & y.number[i]) | (((~x.number[i]) | y.number[i]) & r.number[i])) >> 31;
	}
}

void increment(const BigDouble& x, BigDouble& r) {
	copy(x, r);
	unsigned int carry = 1;
	int i = 0;
	while ((carry == 1) && (i < x.capacity)) {
		r.number[i]++;
		if (r.number[i] != 0)
			carry = 0;
		i++;
	}
}

void mult(const BigDouble& x, const unsigned int y, BigDouble& r) {
	r.nullify();
	unsigned int hh, hl, lh, ll;
	unsigned int temp;
	unsigned int overflow = 0;
	for (int i = 0; i < r.capacity; i++) {
		hh = (x.number[i] >> 16) * (y >> 16);
		hl = (x.number[i] >> 16) * (y & 0xFFFF);
		lh = (x.number[i] & 0xFFFF) * (y >> 16);
		ll = (x.number[i] & 0xFFFF) * (y & 0xFFFF);
		r.number[i] = ll + ((hl + lh) << 16) + overflow;
		temp = overflow;
		overflow = hh + (hl >> 16) + (lh >> 16);
		//for another part of cycle use unsigned int variable "hh" as a buffer for variable "overflow"
		hh = temp;
		temp = (hl & 0xFFFF) + (lh & 0xFFFF);
		overflow = overflow + (temp >> 16);
		temp = (hl + lh) << 16;
		overflow = overflow + (((ll & temp) | ((ll | temp) & (~(ll + temp)))) >> 31);
		temp = temp + ll;
		overflow = overflow + (((hh & temp) | ((hh | temp) & (~(hh + temp)))) >> 31);
	}
}

void mult(const BigDouble& x, const BigDouble& y, BigDouble& r) {
	r.nullify();
	BigDouble temp(0);
	BigDouble aux(0);
	for (int i = 0; i < y.capacity; i++) {
		mult(x, y.number[i], temp);
		leftBlocksShift(temp, i, aux);
		copy(r, temp);
		add(temp, aux, r);
	}
}

void div(const BigDouble& x, const BigDouble& y, BigDouble& q, BigDouble& r) {
	BigDouble temp(0);
	BigDouble aux(0);
	int k = highestNonZeroBit(y) + 1;
	copy(x, r);
	q.nullify();
	while (cmp(r, y) != -1) {
		int t = r.highestNonZeroBit() + 1;
		leftShift(y, t - k, temp);
		if (cmp(r, temp) == -1) {
			t--;
			leftShift(y, t - k, temp);
		}
		copy(r, aux);
		sub(aux, temp, r);
		q.setBit(t - k);
	}
}

void mod(BigDouble& x, BigDouble& n, BigDouble& r) {
	BigDouble temp(0);
	BigDouble aux(0);
	int k = highestNonZeroBit(n) + 1;
	copy(x, r);
	while (cmp(r, n) != -1) {
		int t = r.highestNonZeroBit() + 1;
		leftShift(n, t - k, temp);
		if (cmp(r, temp) == -1) {
			t--;
			leftShift(n, t - k, temp);
		}
		copy(r, aux);
		sub(aux, temp, r);
	}
}

int cmp(const BigDouble &x, const BigDouble& y) {
	int i = x.capacity - 1;
	int res;
	while ((x.number[i] == y.number[i]) && i > -1)
		i--;
	if (i == -1)
		res = 0;
	else if (x.number[i] > y.number[i])
		res = 1;
	else
		res = -1;
	return res;
}

void leftShift(const BigDouble& x, int n, BigDouble& r) {
	r.nullify();
	if (n >= (x.capacity * x.base))
		return;
	int sub_block = n % x.base;
	int blocks = n / x.base;
	int left_lim = x.capacity - blocks;
	int param;
	if (sub_block == 0)
		param = 0;
	else
		param = 0xFFFFFFFF;
	for (int i = 1; i < left_lim; i++) {
		r.number[x.capacity - i] = (x.number[left_lim - i] << sub_block) +
			(param & (x.number[left_lim - (i + 1)] >> (x.base - sub_block)));
	}
	r.number[blocks] = x.number[0] << sub_block;
}

void leftBlocksShift(const BigDouble& x, int n, BigDouble& r) {
	r.nullify();
	if (n >= x.capacity)
		return;
	for (int i = x.capacity - 1; i >= n; i--)
		r.number[i] = x.number[i - n];
}

void rightShift(const BigDouble& x, int n, BigDouble& r) {
	r.nullify();
	if (n >= (x.base * x.capacity))
		return;
	int sub_block = n % x.base;
	int blocks = n / x.base;
	int right_lim = x.capacity - blocks;
	unsigned int param;
	if (sub_block == 0)
		param = 0;
	else
		param = 0xFFFFFFFF >> (x.base - sub_block);
	for (int i = 0; i < right_lim - 1; i++)
		r.number[i] = (x.number[blocks + i] >> sub_block) + ((param & (x.number[blocks + i + 1])) << (x.base - sub_block));
	r.number[right_lim - 1] = x.number[x.capacity - 1] >> sub_block;
}

void rightBlocksShift(const BigDouble& x, int n, BigDouble& r) {
	r.nullify();
	if (n >= x.capacity)
		return;
	for (int i = 0; i < x.capacity - n; i++)
		r.number[i] = x.number[i + n];
}