#include"Interface.h"

int PolynomialPB::irreduciblePolynomial[DENSITY] = { 293, 11, 6, 1, 0 };

PolynomialPB::PolynomialPB() {
	vector = new unsigned int[VECTOR_LENGTH];
	for (int i = 0; i < length; i++)
		vector[i] = 0;
}

PolynomialPB::PolynomialPB(string n) {
	vector = new unsigned int[VECTOR_LENGTH];
	for (int i = 0; i < length; i++)
		vector[i] = 0;
	int j = 0;
	for (int i = n.length() - 1; i > -1; i--) {
		if (n[i] == '1')
			vector[j / base] += (0x1 << (j % base));
		j++;
	}
}

PolynomialPB::~PolynomialPB() {
	delete[] vector;
}

void PolynomialPB::addIdentity() {
	for (int i = 0; i < length; i++)
		vector[i] = 0;
}

void PolynomialPB::multIdentity() {
	vector[0] = 0x1;
	for (int i = 1; i < length; i++)
		vector[i] = 0;
}

void PolynomialPB::set_bin(string n) {
	for (int i = 0; i < length; i++)
		vector[i] = 0;
	int j = 0;
	for (int i = n.length() - 1; i > -1; i--) {
		if (n[i] == '1')
			vector[j / base] += (0x1 << (j % base));
		j++;
	}
}

void PolynomialPB::set_hex(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	for (int i = 0; i < length; i++) {
		vector[i] = 0;
	}
	int index = n.length() - 1;
	int shift = 0;
	int block = 0;
	while (index > -1) {
		int hexadecimal = 0;
		while (codes[hexadecimal] != n[index])
			hexadecimal++;
		vector[block] += hexadecimal << (4 * shift);
		shift++;
		block += (shift / 8);
		shift = shift % 8;
		index--;
	}
}

string PolynomialPB::toString() {
	string result = "";
	for (int i = length - 1; i > -1; i--)
		for (int j = base - 1; j > -1; j--)
			if (((vector[i] >> j) & 0x1) == 1)
				result += '1';
			else
				result += '0';
	result = result.erase(0, (base * length) - fieldExtension);
	return result;
}

string PolynomialPB::toStringHex() {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	string result = "";
	for (int i = length - 1; i > -1; i--)
		for (int j = base / 4 - 1; j > -1; j--)
			result += codes[(vector[i] >> (4 * j)) & 0xF];
	result = result.erase(0, ((base * length) - fieldExtension) / 4);
	return result;
}

int PolynomialPB::degree() {
	int i, j;
	for (i = length - 1; i > -1; i--)
		if (vector[i] != 0)
			break;
	for (j = base - 1; j > -1; j--)
		if (((vector[i] >> j) & 1) != 0)
			break;
	return (base * i + j);
}

void copy(const PolynomialPB& scr, PolynomialPB& dst) {
	for (int i = 0; i < dst.length; i++)
		dst.vector[i] = scr.vector[i];
}

void add(const PolynomialPB& x, const PolynomialPB& y, PolynomialPB& r) {
	for (int i = 0; i < r.length; i++)
		r.vector[i] = x.vector[i] ^ y.vector[i];
}

void mult(const PolynomialPB& x, const PolynomialPB& y, PolynomialPB& r) {
	unsigned int param;
	unsigned int res[2 * VECTOR_LENGTH];
	int DOUBLE_VECTOR_LENGTH = 2 * VECTOR_LENGTH;
	for (int i = 0; i < DOUBLE_VECTOR_LENGTH; i++)
		res[i] = 0;
	for (int i = 0; i < y.length; i++) {
		for (int j = 0; j < y.base; j++) {
			if (((y.vector[i] >> j) & 0x1) == 1) {
				if (j == 0)
					param = 0;
				else
					param = 0xFFFFFFFF;
				res[i] = res[i] ^ (x.vector[0] << j);
				for (int k = 1; k < x.length; k++)
					res[i + k] = res[i + k] ^ ((x.vector[k] << j) + (param & (x.vector[k - 1] >> (x.base - j))));
				res[i + x.length] = res[i + x.length] ^ (param & (x.vector[x.length - 1] >> (x.base - j)));
			}
		}
	}
	for (int i = DOUBLE_VECTOR_LENGTH - 1; i >= x.length - 1; i--) {
		for (int j = x.base - 1; (j > -1) && ((i * x.base + j) >= x.fieldExtension); j--) {
			if (((res[i] >> j) & 0x1) == 1) {
				for (int k = 0; k < DENSITY; k++) {
					int index = (i * x.base + j - x.irreduciblePolynomial[0] + x.irreduciblePolynomial[k]) / x.base;
					int sub_index = (i * x.base + j - x.irreduciblePolynomial[0] + x.irreduciblePolynomial[k]) % x.base;
					res[index] = res[index] ^ (0x1 << sub_index);
				}
			}
		}
	}
	for (int i = 0; i < r.length; i++)
		r.vector[i] = res[i];
}

int trace(const PolynomialPB& x) {
	PolynomialPB temp, res, aux;
	copy(x, res);
	copy(x, temp);
	for (int i = 1; i < x.fieldExtension; i++) {
		sqr(temp, aux);
		copy(res, temp);
		add(temp, aux, res);
		copy(aux, temp);
	}
	return (res.vector[0] & 0x1);
}

void sqr(const PolynomialPB& x, PolynomialPB& r) {
	unsigned int res[2 * VECTOR_LENGTH];
	int DOUBLE_VECTOR_LENGTH = 2 * VECTOR_LENGTH;
	for (int i = 0; i < DOUBLE_VECTOR_LENGTH; i++)
		res[i] = 0;
	for (int i = 0; i < x.fieldExtension; i++) {
		int index = i * 2;
		res[index / x.base] = res[index / x.base] ^ (((x.vector[i / x.base] >> (i % x.base)) & 0x1) << (index % x.base));
	}
	for (int i = DOUBLE_VECTOR_LENGTH - 1; i >= x.length - 1; i--) {
		for (int j = x.base - 1; (j > -1) && ((i * x.base + j) >= x.fieldExtension); j--) {
			if (((res[i] >> j) & 0x1) == 1) {
				for (int k = 0; k < DENSITY; k++) {
					int index = (i * x.base + j - x.irreduciblePolynomial[0] + x.irreduciblePolynomial[k]) / x.base;
					int sub_index = (i * x.base + j - x.irreduciblePolynomial[0] + x.irreduciblePolynomial[k]) % x.base;
					res[index] = res[index] ^ (0x1 << sub_index);
				}
			}
		}
	}
	for (int i = 0; i < r.length; i++)
		r.vector[i] = res[i];
}

void pow(const PolynomialPB& x, unsigned int e, PolynomialPB& r) {
	PolynomialPB temp;
	r.multIdentity();
	for (int i = x.base - 1; i > -1; i--) {
		if (((e >> i) & 0x1) == 1) {
			copy(r, temp);
			mult(x, temp, r);
		}
		if (i != 0) {
			copy(r, temp);
			sqr(temp, r);
		}
	}
}

void pow(const PolynomialPB& x, const BigInteger& e, PolynomialPB& r) {
	PolynomialPB temp;
	r.multIdentity();
	int blocks = highestNonZeroBlock(e);
	for (int i = blocks; i > -1; i--) {
		for (int j = e.base - 1; j > -1; j--) {
			if (((e.number[i] >> j) & 0x1) == 1) {
				copy(r, temp);
				mult(x, temp, r);
			}
			if ((i + j) != 0) {
				copy(r, temp);
				sqr(temp, r);
			}
		}
	}
}

void inverse(const PolynomialPB& x, PolynomialPB& r) {
	BigInteger a(0), b(2), e;
	a.setBit(x.fieldExtension);
	sub(a, b, e);
	pow(x, e, r);
}