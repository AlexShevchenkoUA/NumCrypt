#include"Interface.h"

void binaryAlgorithm(const BigInteger& x, const BigInteger& y, BigInteger& r) {
	BigInteger temp(0);
	BigInteger temp_x(0), temp_y(0);
	BigInteger zero(0);
	int shift = 0;
	copy(x, temp_x);
	copy(y, temp_y);
	while (((temp_x.number[0] & 0x1) == 0) && ((temp_y.number[0] & 0x1) == 0)) {
		rightShift(temp_x, 1, temp);
		copy(temp, temp_x);
		rightShift(temp_y, 1, temp);
		copy(temp, temp_y);
		shift++;
	}
	while ((temp_x.number[0] & 0x1) == 0) {
		copy(temp_x, temp);
		rightShift(temp, 1, temp_x);
	}
	while (cmp(temp_y, zero) != 0) {
		while ((temp_y.number[0] & 0x1) == 0) {
			copy(temp_y, temp);
			rightShift(temp, 1, temp_y);
		}
		if (cmp(temp_x, temp_y) != -1) {
			copy(temp_y, temp);
			sub(temp_x, temp, temp_y);
			copy(temp, temp_x);
		}
		else {
			copy(temp_y, temp);
			sub(temp, temp_x, temp_y);
		}
	}
	leftShift(temp_x, shift, r);
}

void linearCongruencesSystem(const BigInteger* b, const BigInteger* n, int size, BigInteger& r) {
	BigInteger temp, aux;
	BigInteger M(1);
	BigInteger* N = new BigInteger[size];
	BigInteger* m = new BigInteger[size];
	for (int i = 0; i < size; i++) {
		mult(M, n[i], temp);
		copy(temp, M);
	}
	for (int i = 0; i < size; i++) 
		div(M, n[i], m[i], temp);
	for (int i = 0; i < size; i++)
		inverseElement(m[i], n[i], N[i]);
	r.nullify();
	for (int i = 0; i < size; i++) {
		modMult(b[i], N[i], M, temp);
		modMult(temp, m[i], M, aux);
		copy(r, temp);
		modAdd(temp, aux, M, r);
	}
	delete[] m;
	delete[] N;
}

void random(int blocks, BigInteger& r) {
	unsigned int temp;
	r.nullify();
	srand((unsigned int)time(NULL));
	for (int i = 0; i < blocks; i++) {
		temp = 0;
		for (int j = 0; j < 4; j++) {
			unsigned int aux = rand() % 0x100;
			temp += (aux << (8 * j));
		}
		r.number[i] = temp;
	}
}

void MillerRabinPrimalityTest(const BigInteger& n, int& r) {
	int rounds = MILLER_RABIN_PRIMALITY_TEST_CONST;
	int s = 0;
	BigInteger d, a, x;
	BigInteger temp;
	BigInteger identity(1), nidentity;
	sub(n, identity, d);
	copy(d, nidentity);
	int random_blocks = nidentity.highestNonZeroBlock() + 1;
	while (((d.number[s / d.base] >> (s % d.base)) & 0x1) == 0)
		s++;
	rightShift(nidentity, s, d);
	for (int i = 0; i < rounds; i++) {
		random(random_blocks, a);
		modPowBarrett(a, d, n, x);
		if ((cmp(x, identity) == 0) || (cmp(x, nidentity) == 0))
			continue;
		bool flag = false;
		for (int j = 0; j < s - 1; j++) {
			copy(x, temp);
			modSqr(temp, n, x);
			if (cmp(x, identity) == 0) {
				r = 0;
				return;
			}
			if (cmp(x, nidentity) == 0) {
				flag = true;
				break;
			}
		}
		if (flag)
			continue;
		r = 0;
		return;
	}
	r = 1;
}