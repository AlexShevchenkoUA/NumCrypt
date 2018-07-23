#include"Interface.h"

bool PolynomialNB::status = false;

unsigned int PolynomialNB::multiplicativeMatrix[FIELD_EXTENSION][2][2] = { 0 };

void PolynomialNB::calculatingMultiplicativeMatrix() {
	BigInteger p(2 * FIELD_EXTENSION + 1);
	BigInteger temp(0), aux(0);
	BigInteger identity(1), nidentity(2 * FIELD_EXTENSION);
	BigInteger res(0);
	BigInteger cache[FIELD_EXTENSION];
	for (int i = 0; i < FIELD_EXTENSION; i++) {
		cache[i].nullify();
		temp.nullify();
		temp.setBit(i);
		mod(temp, p, cache[i]);
	}
	multiplicativeMatrix[0][0][0] = 0;
	multiplicativeMatrix[0][0][1] = 1;
	multiplicativeMatrix[0][1][0] = 0;
	multiplicativeMatrix[0][1][1] = 1;
	int flag;
	for (int i = 1; i < FIELD_EXTENSION; i++) {
		flag = 0;
		for (int j = 0; j < FIELD_EXTENSION; j++) {
			modAdd(cache[i], cache[j], p, res);
			if ((cmp(identity, res) * cmp(nidentity, res)) == 0) {
				multiplicativeMatrix[i][flag][0] = j / base;
				multiplicativeMatrix[i][flag][1] = j % base;
				flag++;
				if (flag == 2)
					break;
				continue;
 			}
			modSub(cache[i], cache[j], p, res);
			if ((cmp(identity, res) * cmp(nidentity, res)) == 0) {
				multiplicativeMatrix[i][flag][0] = j / base;
				multiplicativeMatrix[i][flag][1] = j % base;
				flag++;
				if (flag == 2)
					break;
				continue;
			}
		}
	}
}

PolynomialNB::PolynomialNB() {
	if (!status) {
		calculatingMultiplicativeMatrix();
		status = true;
	}
	vector = new unsigned int[VECTOR_LENGTH];
	for (int i = 0; i < length; i++)
		vector[i] = 0;
}

PolynomialNB::PolynomialNB(string n) {
	if (!status) {
		calculatingMultiplicativeMatrix();
		status = true;
	}
	vector = new unsigned int[VECTOR_LENGTH];
	for (int i = 0; i < length; i++)
		vector[i] = 0;
	for (int i = 0; i < n.length(); i++) 
		if (n[i] == '1')
			vector[i / base] += (0x1 << (i % base));
}

PolynomialNB::~PolynomialNB() {
	delete[] vector;
}

void PolynomialNB::addIdentity() {
	for (int i = 0; i < length; i++)
		vector[i] = 0;
}

void PolynomialNB::multIdentity() {
	for (int i = 0; i < length; i++)
		vector[i] = 0xFFFFFFFF;
	vector[length - 1] = vector[length - 1] >> (base * length - fieldExtension);
}

void PolynomialNB::set_bin(string n) {
	for (int i = 0; i < length; i++)
		vector[i] = 0;
	for (int i = 0; i < n.length(); i++) {
		if (n[i] == '1')
			vector[i / base] += (0x1 << (i % base));
	}
}

void PolynomialNB::set_hex(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	string binary = "";
	for (int i = 0; i < n.length(); i++) {
		unsigned int hex = 0;
		while (codes[hex] != n[i])
			hex++;
		for (int j = 3; j > -1; j--) {
			if (((hex >> j) & 0x1) == 1)
				binary += '1';
			else
				binary += '0';
		}
	}
	binary.erase(0, 4 - (fieldExtension % 4));
	for (int i = 0; i < length; i++)
		vector[i] = 0;
	for (int i = 0; i < binary.length(); i++) {
		if (binary[i] == '1')
			vector[i / base] += (0x1 << (i % base));
	}
}

string PolynomialNB::toString() {
	string result = "";
	for (int i = 0; i < length; i++)
		for (int j = 0; j < base; j++)
			if (((vector[i] >> j) & 0x1) == 1)
				result += '1';
			else
				result += '0';
	result = result.erase(fieldExtension, (base * length) - fieldExtension);
	return result;
}

string PolynomialNB::toStringHex() {
	BigInteger res, temp, aux;
	for (int i = 0; i < length; i++) {
		unsigned int t = vector[i];
		//Reversing bits
		t = ((t & 0x55555555) << 1) | ((t & 0xAAAAAAAA) >> 1);
		t = ((t & 0x33333333) << 2) | ((t & 0xCCCCCCCC) >> 2);
		t = ((t & 0x0F0F0F0F) << 4) | ((t & 0xF0F0F0F0) >> 4);
		t = ((t & 0x00FF00FF) << 8) | ((t & 0xFF00FF00) >> 8);
		t = ((t & 0x0000FFFF) << 16) | ((t & 0xFFFF0000) >> 16);
		leftBlocksShift(res, 1, temp);
		aux.set(t);
		add(temp, aux, res);
	}
	rightShift(res, length * base - fieldExtension, temp);
	string result = temp.toString();
	while (result.length() < ((fieldExtension / 4) + 1))
		result = "0" + result;
	return result; 
}

void copy(const PolynomialNB& scr, PolynomialNB& dst) {
	for (int i = 0; i < dst.length; i++)
		dst.vector[i] = scr.vector[i];
}

void add(const PolynomialNB& x, const PolynomialNB& y, PolynomialNB& r) {
	for (int i = 0; i < r.length; i++)
		r.vector[i] = x.vector[i] ^ y.vector[i];
}

void mult(const PolynomialNB& x, const PolynomialNB& y, PolynomialNB& r) {
	unsigned int aux = 0;
	r.addIdentity();
	PolynomialNB temp_x, temp_y;
	for (int i = 0; i < r.fieldExtension; i++) {
		leftRotation(x, i, temp_x);
		leftRotation(y, i, temp_y);
		tensorMult(temp_x, temp_y, aux);
		r.vector[i / r.base] += (aux << (i % r.base));
	}
}

void tensorMult(const PolynomialNB& x, const PolynomialNB& y, unsigned int& r) {
	unsigned int temp;
	r = ((x.vector[x.multiplicativeMatrix[0][0][0]] >> (x.multiplicativeMatrix[0][0][1])) & 0x1) & (y.vector[0] & 0x1);
	for (int i = 1; i < x.fieldExtension; i++) {
		temp = ((x.vector[x.multiplicativeMatrix[i][0][0]] >> (x.multiplicativeMatrix[i][0][1])) & 0x1) ^
			((x.vector[x.multiplicativeMatrix[i][1][0]] >> (x.multiplicativeMatrix[i][1][1])) & 0x1);
		temp = temp & ((y.vector[i / y.base] >> (i % y.base)) & 0x1);
		r = r ^ temp;
	}
}

int trace(PolynomialNB& x) {
	unsigned int res = 0;
	for (int i = 0; i < x.length; i++) {
		int temp = x.vector[i] ^ (x.vector[i] >> 16);
		temp = temp ^ (temp >> 8);
		temp = temp ^ (temp >> 4);
		temp = temp ^ (temp >> 2);
		temp = temp ^ (temp >> 1);
		res = res ^ (temp & 0x1);
	}
	return res;
}

void sqr(const PolynomialNB& x, PolynomialNB& r) {
	r.addIdentity();
	r.vector[0] = (x.vector[0] << 1) | (x.vector[x.length - 1] >> ((x.fieldExtension % x.base) + 1));
	r.vector[r.length - 1] = ((x.vector[x.length - 1] << 1) | (x.vector[x.length - 2] >> (x.base - 1))) & 
								(0xFFFFFFFF >> (x.base * x.length - x.fieldExtension));
	for (int i = 1; i < x.length - 1; i++)
		r.vector[i] = (x.vector[i] << 1) | (x.vector[i - 1] >> (x.base - 1));
}

void pow(const PolynomialNB& x, unsigned int e, PolynomialNB& r) {
	PolynomialNB temp;
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

void pow(const PolynomialNB& x, const BigInteger& e, PolynomialNB& r) {
	PolynomialNB temp;
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

void inverse(const PolynomialNB& x, PolynomialNB& r) {
	PolynomialNB b, temp, aux;
	r.addIdentity();
	unsigned int e = x.fieldExtension - 1;
	int k = 1;
	int t = x.base - 1;
	while (((e >> t) & 0x1) == 0) {
		t--;
	}
	copy(x, b);
	for (int i = t - 1; i > -1; i--) {
		copy(b, aux);
		rightRotation(b, k, temp);
		mult(temp, aux, b);
		k = 2 * k;
		if (((e >> i) & 0x1) == 0x1) {
			rightRotation(b, 1, temp);
			mult(x, temp, b);
			k = k + 1;
		}
	}
	rightRotation(b, 1, r); 
}

void leftRotation(const PolynomialNB& x, int n, PolynomialNB& r) {
	r.addIdentity();
	rightRotation(x, x.fieldExtension - (n % x.fieldExtension), r);
}

void rightRotation(const PolynomialNB& x, int n, PolynomialNB& r) {
	r.addIdentity();
	int rotation = n % x.fieldExtension;
	unsigned int param;
	unsigned int res[2 * VECTOR_LENGTH];
	unsigned int temp[VECTOR_LENGTH]; 
	int DOUBLE_VECTOR_LENGTH = 2 * VECTOR_LENGTH;
	for (int i = 0; i < DOUBLE_VECTOR_LENGTH; i++)
		res[i] = 0;
	for (int i = 0; i < VECTOR_LENGTH; i++)
		temp[i] = 0;
	int blocks = rotation / x.base;
	int shift = rotation % x.base;
	if (shift == 0)
		param = 0;
	else
		param = 0xFFFFFFFF;
	res[blocks] = res[blocks] ^ (x.vector[0] << shift);
	for (int k = 1; k < x.length; k++)
		res[blocks + k] = (x.vector[k] << shift) + (param & (x.vector[k - 1] >> (x.base - shift)));
	res[blocks + x.length] = param & (x.vector[x.length - 1] >> (x.base - shift));
	int j = 0;
	for (int i = x.length - 1 - j; j < VECTOR_LENGTH; i++) {
		temp[j] = res[i];
		j++;
	}
	param = 0xFFFFFFFF >> (x.length * x.base - x.fieldExtension);
	int sub_block = x.fieldExtension - (x.length - 1) * x.base;
	for (int i = 0; i < VECTOR_LENGTH - 1; i++) {
		unsigned int aux = (temp[i] >> sub_block) | ((temp[i + 1] & param) << (x.base - sub_block));
		res[i] = res[i] | aux;
	}
	res[VECTOR_LENGTH - 1] = res[VECTOR_LENGTH - 1] | (temp[VECTOR_LENGTH - 1] >> sub_block);
	for (int i = 0; i < VECTOR_LENGTH; i++) {
		r.vector[i] = res[i];
	}
}

/*unsigned int PolynomialNB::multiplicativeMatrix[FIELD_EXTENSION][2][2] = { 
	{ { 0, 1}, { 0, 1} },
	{ { 0, 0}, { 5, 25} },
	{ { 5, 25}, { 6, 13} },
	{ { 1, 21}, { 2, 13} },
	{ { 1, 6}, { 3, 1} },
	{ { 1, 29}, { 8, 12} },
	{ { 4, 2}, { 6, 25} },
	{ { 7, 4}, { 8, 21} },
	{ { 4, 7}, { 5, 21} },
	{ { 6, 10}, { 7, 11} },
	{ { 1, 4}, { 2, 11} },
	{ { 2, 18}, { 3, 0} },
	{ { 1, 22}, { 5, 27} },
	{ { 3, 21}, { 8, 6} },
	{ { 6, 20}, { 6, 31} },
	{ { 4, 22}, { 5, 8} },
	{ { 0, 23}, { 4, 19} },
	{ { 1, 13}, { 7, 20} },
	{ { 4, 16}, { 4, 31} },
	{ { 3, 11}, { 5, 30} },
	{ { 3, 13}, { 3, 15} },
	{ { 3, 24}, { 5, 1} },
	{ { 4, 18}, { 5, 18} },
	{ { 0, 16}, { 4, 23} },
	{ { 7, 17}, { 7, 20} },
	{ { 0, 30}, { 2, 0} },
	{ { 2, 22}, { 8, 27} },
	{ { 5, 4}, { 7, 2} },
	{ { 4, 14}, { 8, 20} },
	{ { 1, 31}, { 4, 3} },
	{ { 0, 25}, { 7, 18} },
	{ { 1, 12}, { 2, 22} },
	{ { 4, 20}, { 5, 10} },
	{ { 2, 28}, { 3, 5} },
	{ { 8, 8}, { 9, 1} },
	{ { 2, 10}, { 3, 18} },
	{ { 0, 10}, { 7, 12} },
	{ { 3, 0}, { 7, 8} },
	{ { 0, 4}, { 2, 14} },
	{ { 8, 2}, { 8, 12} },
	{ { 6, 28}, { 7, 31} },
	{ { 5, 14}, { 6, 3} },
	{ { 7, 26}, { 8, 25} },
	{ { 2, 21}, { 7, 15} },
	{ { 0, 31}, { 7, 19} },
	{ { 0, 17}, { 4, 20} },
	{ { 4, 31}, { 5, 7} },
	{ { 3, 23}, { 7, 0} },
	{ { 3, 16}, { 6, 0} },
	{ { 2, 2}, { 2, 9} },
	{ { 2, 30}, { 9, 2} },
	{ { 2, 17}, { 5, 23} },
	{ { 2, 12}, { 6, 12} },
	{ { 0, 3}, { 5, 26} },
	{ { 0, 12}, { 3, 1} },
	{ { 4, 11}, { 8, 6} },
	{ { 3, 29}, { 9, 0} },
	{ { 2, 29}, { 5, 12} },
	{ { 2, 3}, { 6, 2} },
	{ { 8, 0}, { 8, 4} },
	{ { 4, 10}, { 8, 11} },
	{ { 0, 5}, { 3, 2} },
	{ { 4, 2}, { 4, 29} },
	{ { 0, 29}, { 4, 15} },
	{ { 0, 25}, { 7, 21} },
	{ { 2, 8}, { 8, 27} },
	{ { 1, 17}, { 6, 1} },
	{ { 1, 26}, { 2, 30} },
	{ { 7, 7}, { 8, 4} },
	{ { 3, 20}, { 7, 13} },
	{ { 2, 20}, { 5, 28} },
	{ { 5, 16}, { 8, 26} },
	{ { 2, 1}, { 7, 22} },
	{ { 1, 17}, { 3, 17} },
	{ { 1, 3}, { 9, 2} },
	{ { 0, 10}, { 6, 11} },
	{ { 1, 20}, { 2, 18} },
	{ { 0, 3}, { 6, 14} },
	{ { 1, 6}, { 7, 9} },
	{ { 6, 9}, { 8, 2} },
	{ { 5, 22}, { 7, 6} },
	{ { 1, 19}, { 2, 31} },
	{ { 0, 11}, { 2, 12} },
	{ { 5, 27}, { 6, 14} },
	{ { 2, 6}, { 7, 14} },
	{ { 1, 11}, { 8, 26} },
	{ { 0, 26}, { 0, 31} },
	{ { 5, 4}, { 5, 10} },
	{ { 2, 26}, { 8, 18} },
	{ { 4, 13}, { 8, 17} },
	{ { 2, 24}, { 5, 5} },
	{ { 3, 4}, { 8, 17} },
	{ { 1, 1}, { 5, 11} },
	{ { 1, 25}, { 9, 1} },
	{ { 1, 18}, { 2, 3} },
	{ { 2, 17}, { 7, 7} },
	{ { 0, 11}, { 1, 5} },
	{ { 0, 4}, { 1, 22} },
	{ { 1, 29}, { 4, 11} },
	{ { 4, 29}, { 8, 16} },
	{ { 2, 27}, { 5, 6} },
	{ { 1, 1}, { 4, 21} },
	{ { 6, 21}, { 8, 8} },
	{ { 3, 26}, { 7, 24} },
	{ { 6, 18}, { 8, 24} },
	{ { 5, 15}, { 6, 30} },
	{ { 3, 22}, { 5, 29} },
	{ { 0, 19}, { 5, 0} },
	{ { 3, 13}, { 3, 14} },
	{ { 0, 20}, { 3, 12} },
	{ { 3, 12}, { 5, 1} },
	{ { 0, 20}, { 5, 31} },
	{ { 1, 16}, { 3, 24} },
	{ { 2, 9}, { 7, 23} },
	{ { 1, 3}, { 8, 9} },
	{ { 4, 9}, { 7, 12} },
	{ { 2, 5}, { 8, 5} },
	{ { 0, 13}, { 5, 28} },
	{ { 3, 10}, { 6, 31} },
	{ { 1, 15}, { 5, 0} },
	{ { 0, 21}, { 3, 16} },
	{ { 5, 18}, { 7, 23} },
	{ { 3, 7}, { 6, 22} },
	{ { 4, 27}, { 6, 18} },
	{ { 8, 15}, { 8, 31} },
	{ { 1, 24}, { 4, 12} },
	{ { 5, 12}, { 8, 19} },
	{ { 7, 3}, { 8, 29} },
	{ { 6, 26}, { 8, 14} },
	{ { 4, 28}, { 6, 24} },
	{ { 0, 6}, { 1, 30} },
	{ { 0, 29}, { 8, 21} },
	{ { 4, 25}, { 7, 18} },
	{ { 7, 16}, { 7, 28} },
	{ { 4, 24}, { 5, 20} },
	{ { 0, 8}, { 8, 22} },
	{ { 6, 6}, { 7, 11} },
	{ { 3, 19}, { 8, 10} },
	{ { 1, 28}, { 8, 5} },
	{ { 1, 23}, { 3, 2} },
	{ { 3, 29}, { 8, 16} },
	{ { 2, 25}, { 8, 19} },
	{ { 0, 28}, { 5, 5} },
	{ { 1, 31}, { 4, 30} },
	{ { 0, 18}, { 7, 21} },
	{ { 5, 17}, { 5, 30} },
	{ { 0, 22}, { 5, 2} },
	{ { 0, 16}, { 5, 9} },
	{ { 1, 0}, { 1, 13} },
	{ { 3, 5}, { 5, 7} },
	{ { 0, 15}, { 6, 21} },
	{ { 0, 23}, { 5, 19} },
	{ { 4, 6}, { 7, 17} },
	{ { 4, 4}, { 8, 22} },
	{ { 6, 17}, { 7, 28} },
	{ { 3, 27}, { 6, 23} },
	{ { 4, 1}, { 8, 15} },
	{ { 1, 30}, { 3, 3} },
	{ { 4, 15}, { 5, 6} },
	{ { 0, 18}, { 1, 14} },
	{ { 3, 11}, { 3, 23} },
	{ { 0, 21}, { 3, 14} },
	{ { 4, 18}, { 5, 31} },
	{ { 5, 9}, { 7, 1} },
	{ { 0, 27}, { 2, 23} },
	{ { 2, 26}, { 4, 14} },
	{ { 3, 4}, { 4, 30} },
	{ { 1, 14}, { 4, 21} },
	{ { 0, 15}, { 7, 0} },
	{ { 4, 19}, { 5, 3} },
	{ { 1, 0}, { 2, 23} },
	{ { 2, 28}, { 8, 18} },
	{ { 1, 25}, { 3, 30} },
	{ { 6, 2}, { 8, 29} },
	{ { 1, 9}, { 6, 29} },
	{ { 3, 9}, { 8, 25} },
	{ { 2, 7}, { 5, 29} },
	{ { 4, 17}, { 7, 22} },
	{ { 0, 22}, { 3, 25} },
	{ { 4, 23}, { 6, 22} },
	{ { 4, 6}, { 7, 29} },
	{ { 0, 8}, { 7, 5} },
	{ { 2, 16}, { 6, 10} },
	{ { 1, 19}, { 9, 3} },
	{ { 6, 12}, { 9, 4} },
	{ { 0, 1}, { 0, 2} },
	{ { 1, 21}, { 6, 13} },
	{ { 0, 12}, { 2, 19} },
	{ { 2, 6}, { 3, 21} },
	{ { 3, 10}, { 5, 16} },
	{ { 0, 19}, { 4, 17} },
	{ { 3, 15}, { 5, 2} },
	{ { 1, 16}, { 7, 1} },
	{ { 2, 2}, { 8, 28} },
	{ { 1, 26}, { 5, 13} },
	{ { 1, 9}, { 8, 0} },
	{ { 6, 8}, { 7, 26} },
	{ { 6, 16}, { 7, 10} },
	{ { 4, 8}, { 8, 23} },
	{ { 7, 25}, { 8, 10} },
	{ { 6, 4}, { 8, 1} },
	{ { 2, 15}, { 7, 10} },
	{ { 0, 9}, { 5, 22} },
	{ { 2, 11}, { 9, 3} },
	{ { 1, 20}, { 5, 24} },
	{ { 0, 2}, { 5, 26} },
	{ { 2, 13}, { 2, 19} },
	{ { 7, 9}, { 7, 14} },
	{ { 6, 5}, { 7, 27} },
	{ { 4, 26}, { 8, 23} },
	{ { 3, 8}, { 3, 27} },
	{ { 6, 30}, { 8, 31} },
	{ { 0, 14}, { 8, 7} },
	{ { 3, 6}, { 4, 22} },
	{ { 3, 26}, { 5, 19} },
	{ { 4, 27}, { 7, 29} },
	{ { 4, 1}, { 6, 27} },
	{ { 0, 6}, { 8, 13} },
	{ { 4, 0}, { 7, 4} },
	{ { 6, 24}, { 7, 30} },
	{ { 1, 8}, { 8, 13} },
	{ { 5, 14}, { 8, 30} },
	{ { 3, 9}, { 6, 19} },
	{ { 0, 14}, { 3, 22} },
	{ { 1, 15}, { 5, 8} },
	{ { 5, 3}, { 6, 0} },
	{ { 0, 27}, { 8, 28} },
	{ { 3, 31}, { 8, 20} },
	{ { 0, 7}, { 6, 26} },
	{ { 5, 21}, { 7, 30} },
	{ { 2, 16}, { 8, 3} },
	{ { 2, 4}, { 2, 31} },
	{ { 1, 5}, { 7, 13} },
	{ { 2, 14}, { 6, 15} },
	{ { 6, 5}, { 6, 9} },
	{ { 0, 9}, { 4, 8} },
	{ { 1, 4}, { 3, 19} },
	{ { 2, 5}, { 7, 8} },
	{ { 2, 20}, { 6, 15} },
	{ { 1, 11}, { 7, 27} },
	{ { 4, 5}, { 7, 19} },
	{ { 0, 24}, { 4, 24} },
	{ { 0, 30}, { 4, 4} },
	{ { 1, 12}, { 7, 16} },
	{ { 0, 17}, { 0, 24} },
	{ { 2, 0}, { 4, 16} },
	{ { 2, 8}, { 5, 17} },
	{ { 3, 17}, { 3, 25} },
	{ { 3, 7}, { 8, 9} },
	{ { 6, 7}, { 8, 24} },
	{ { 1, 10}, { 6, 4} },
	{ { 6, 16}, { 7, 15} },
	{ { 4, 5}, { 4, 26} },
	{ { 5, 20}, { 6, 23} },
	{ { 6, 27}, { 7, 5} },
	{ { 1, 8}, { 8, 3} },
	{ { 1, 27}, { 6, 3} },
	{ { 6, 8}, { 8, 11} },
	{ { 1, 7}, { 2, 15} },
	{ { 7, 6}, { 7, 31} },
	{ { 1, 27}, { 2, 4} },
	{ { 3, 20}, { 4, 10} },
	{ { 0, 13}, { 1, 23} },
	{ { 6, 20}, { 9, 0} },
	{ { 1, 2}, { 3, 6} },
	{ { 3, 18}, { 7, 24} },
	{ { 4, 9}, { 6, 7} },
	{ { 1, 28}, { 8, 1} },
	{ { 0, 5}, { 1, 7} },
	{ { 6, 25}, { 6, 28} },
	{ { 4, 0}, { 8, 30} },
	{ { 3, 28}, { 4, 28} },
	{ { 3, 3}, { 4, 12} },
	{ { 2, 25}, { 2, 27} },
	{ { 2, 24}, { 5, 11} },
	{ { 3, 30}, { 4, 13} },
	{ { 0, 28}, { 7, 3} },
	{ { 0, 7}, { 4, 3} },
	{ { 4, 7}, { 4, 25} },
	{ { 6, 6}, { 6, 17} },
	{ { 3, 8}, { 7, 25} },
	{ { 1, 10}, { 5, 15} },
	{ { 2, 7}, { 2, 21} },
	{ { 0, 26}, { 2, 1} },
	{ { 6, 1}, { 7, 2} },
	{ { 3, 31}, { 5, 13} },
	{ { 6, 29}, { 8, 14} },
	{ { 3, 28}, { 6, 19} },
	{ { 1, 24}, { 8, 7} },
	{ { 1, 2}, { 2, 29} },
	{ { 1, 18}, { 2, 10} },
	{ { 5, 23}, { 6, 11} },
	{ { 5, 24}, { 9, 4} }
	};
*/
