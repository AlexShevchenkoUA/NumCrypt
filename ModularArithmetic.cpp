#include"Interface.h"

void gcd(const BigInteger& x, const BigInteger& y, BigInteger& r) {
	BigInteger remainder(0);
	BigInteger temp_1(0);
	BigInteger temp_2(0);
	BigInteger zero(0);
	copy(x, temp_1);
	copy(y, temp_2);
	while (true) {
		mod(temp_1, temp_2, remainder);
		if (cmp(remainder, zero) == 0)
			break;
		copy(temp_2, temp_1);
		copy(remainder, temp_2);
	}
	copy(temp_2, r);
}

void lcp(const BigInteger& x, const BigInteger& y, BigDouble& r) {
	BigInteger itemp;
	BigDouble dtemp_1;
	BigDouble dtemp_2;
	BigDouble dtemp_3;
	gcd(x, y, itemp);
	mult(x, y, dtemp_1);
	copy(itemp, dtemp_2);
	div(dtemp_1, dtemp_2, r, dtemp_3);
}

void inverseElement(const BigInteger& x, const BigInteger& n, BigInteger& r) {
	BigInteger coef_1(1);
	BigInteger coef_2(0);
	BigInteger a, b, q, rem, temp;
	temp.set(1);
	if (cmp(x, temp) == 0) {
		r.set(1);
		return;
	}
	q.nullify();
	rem.nullify();
	r.nullify();
	copy(n, a);
	copy(x, b);
	int i = 0;
	while (true) {
		i++;
		div(a, b, q, rem);
		mult(coef_1, q, temp);
		add(coef_2, temp, r);
		copy(coef_1, coef_2);
		copy(r, coef_1);
		copy(b, a);
		copy(rem, b);
		temp.set(1);
		if (cmp(rem, temp) == 0) //verify
			break;
	}
	if ((i % 2) == 1) {
		copy(r, temp);
		sub(n, temp, r);
	}
}

void modAdd(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r) {
	r.nullify();
	BigDouble temp_1(0), temp_2(0), temp_3(0);
	add(x, y, temp_1);
	copy(n, temp_2);
	mod(temp_1, temp_2, temp_3);
	copy(temp_3, r);
}

void modSub(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r) {
	r.nullify();
	BigInteger temp_1(0);
	BigInteger temp_2(0);
	mod(x, n, temp_1);
	mod(y, n, temp_2);
	if (cmp(temp_1, temp_2) != -1)
		sub(temp_1, temp_2, r);
	else {
		sub(temp_2, temp_1, r);
		copy(r, temp_1);
		sub(n, temp_1, r);
	}
}

void modMult(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r) {
	BigDouble temp_1, temp_2, temp_3;
	mult(x, y, temp_1);
	copy(n, temp_2);
	mod(temp_1, temp_2, temp_3);
	copy(temp_3, r);
}

void modSqr(const BigInteger& x, const BigInteger& n, BigInteger& r) {
	BigDouble temp_1, temp_2, temp_3;
	mult(x, x, temp_1);
	copy(n, temp_2);
	mod(temp_1, temp_2, temp_3);
	copy(temp_3, r);
}

void barrettReduction(const BigDouble& x, const BigInteger& n, const BigDouble& m, BigInteger& r) {
	BigDouble res(0);
	BigDouble temp(0);
	BigDouble aux(0);
	BigInteger itemp(0);
	int k = highestNonZeroBlock(n) + 1;
	rightBlocksShift(x, k - 1, temp);
	mult(temp, m, res);
	rightBlocksShift(res, k + 1, aux);
	copy(n, temp);
	mult(aux, temp, res);
	sub(x, res, temp);
	copy(temp, r);
	while (cmp(r, n) != -1) {
		copy(r, itemp);
		sub(itemp, n, r);
	}
}

void barrettConst(const BigInteger& n, BigDouble& m) {
	int k = highestNonZeroBlock(n) + 1;
	BigDouble dtemp_n(0);
	BigDouble temp(0);
	BigDouble aux(0);
	copy(n, dtemp_n);
	if (k == n.capacity) {
		temp.maxNumber();
		div(temp, dtemp_n, m, aux);
		increment(aux, temp);
		if (cmp(temp, dtemp_n) == 0) {
			copy(m, temp);
			increment(temp, m);
		}
	}
	else {
		temp.setBit(2 * k * n.base);
		div(temp, dtemp_n, m, aux);
	}
}

void modPowBarrett(const BigInteger& x, const BigInteger& e, const BigInteger& n, BigInteger& r) {
	r.set(1);
	BigDouble m(0);
	BigInteger itemp(0);
	BigDouble dtemp(0);
	mod(x, n, itemp);
	barrettConst(n, m);
	int blocks = highestNonZeroBlock(e) + 1;
	for (int i = 0; i < blocks; i++) {
		for (int j = 0; j < e.base; j++) {
			if (((e.number[i] >> j) & 0x1) == 1) {
				mult(r, itemp, dtemp);
				barrettReduction(dtemp, n, m, r);
			}
			mult(itemp, itemp, dtemp);
			barrettReduction(dtemp, n, m, itemp);
		}
	}
}