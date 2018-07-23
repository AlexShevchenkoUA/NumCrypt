#include"Interface.h"

//class BigInteger

BigInteger::BigInteger() {
	number = new unsigned int[capacity];
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
}

BigInteger::BigInteger(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	number = new unsigned int[capacity];
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
	int index = n.length() - 1;
	int shift = 0;
	int block = 0;
	while (index > -1) {
		int hexadecimal = 0;
		while (codes[hexadecimal] != n[index])
			hexadecimal++;
		number[block] += hexadecimal << (4 * shift);
		shift++;
		block += (shift / 8);
		shift = shift % 8;
		index--;
	}
}

BigInteger::BigInteger(unsigned int n) {
	number = new unsigned int[capacity];
	number[0] = n;
	for (int i = 1; i < capacity; i++) {
		number[i] = 0;
	}
}

BigInteger::~BigInteger() {
	delete[] number;
}

string BigInteger::toString() {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	string result = "";
	int length = capacity - 1;
	while ((number[length] == 0) && (length > 0))
		length--;
	for (int i = length; i > -1; i--)
		for (int j = base / 4 - 1; j > -1; j--)
			result += codes[(number[i] >> (4 * j)) & 0xF];
	while ((result.length() > 1) && (result[0] == '0'))
		result = result.erase(0, 1);
	return result;
}

void BigInteger::set(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
	int index = n.length() - 1;
	int shift = 0;
	int block = 0;
	while (index > -1) {
		int hexadecimal = 0;
		while (codes[hexadecimal] != n[index])
			hexadecimal++;
		number[block] += hexadecimal << (4 * shift);
		shift++;
		block += (shift / 8);
		shift = shift % 8;
		index--;
	}
}

void BigInteger::set(unsigned int n) {
	number[0] = n;
	for (int i = 1; i < capacity; i++) {
		number[i] = 0;
	}
}
void BigInteger::setBit(int n) {
	int block = n / base;
	unsigned int bit_block = 0x1 << (n % base);
	number[block] = number[block] | bit_block;
}

void BigInteger::nullify() {
	for (int i = 0; i < capacity; i++)
		number[i] = 0;
}

void BigInteger::maxNumber() {
	for (int i = 0; i < capacity; i++)
		number[i] = 0xFFFFFFFF;
}

int BigInteger::highestNonZeroBit() {
	int i, j;
	for (i = capacity - 1; i > 0; i--) {
		if (number[i] != 0)
			break;
	}
	for (j = base - 1; j > 0; j--) {
		if (((number[i] >> j) & 1) != 0)
			break;
	}
	return (base * i + j);
}

int BigInteger::highestNonZeroBlock() {
	int i;
	for (i = capacity - 1; i > 0; i--)
		if (number[i] != 0)
			break;
	return i;
}

//class BigDouble

BigDouble::BigDouble() {
	number = new unsigned int[capacity];
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
}

BigDouble::BigDouble(unsigned int n) {
	number = new unsigned int[capacity];
	number[0] = n;
	for (int i = 1; i < capacity; i++) {
		number[i] = 0;
	}
}

BigDouble::BigDouble(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	number = new unsigned int[capacity];
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
	int index = n.length() - 1;
	int shift = 0;
	int block = 0;
	while (index > -1) {
		int hexadecimal = 0;
		while (codes[hexadecimal] != n[index])
			hexadecimal++;
		number[block] += hexadecimal << (4 * shift);
		shift++;
		block += (shift / 8);
		shift = shift % 8;
		index--;
	}
}

BigDouble::~BigDouble() {
	delete[] number;
}

string BigDouble::toString() {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	string result = "";
	int length = capacity - 1;
	while ((number[length] == 0) && (length > 0))
		length--;
	for (int i = length; i > -1; i--)
		for (int j = base / 4 - 1; j > -1; j--)
			result += codes[(number[i] >> (4 * j)) & 0xF];
	while ((result.length() > 1) && (result[0] == '0'))
		result = result.erase(0, 1);
	return result;
}

void BigDouble::set(unsigned int n) {
	number[0] = n;
	for (int i = 1; i < capacity; i++) {
		number[i] = 0;
	}
}

void BigDouble::set(string n) {
	char codes[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
	for (int i = 0; i < capacity; i++) {
		number[i] = 0;
	}
	int index = n.length() - 1;
	int shift = 0;
	int block = 0;
	while (index > -1) {
		int hexadecimal = 0;
		while (codes[hexadecimal] != n[index])
			hexadecimal++;
		number[block] += hexadecimal << (4 * shift);
		shift++;
		block += (shift / 8);
		shift = shift % 8;
		index--;
	}
}

void BigDouble::setBit(int n) {
	int block = n / base;
	unsigned int bit_block = 0x1 << (n % base);
	number[block] = number[block] | bit_block;
}

void BigDouble::nullify() {
	for (int i = 0; i < capacity; i++)
		number[i] = 0;
}

void BigDouble::maxNumber() {
	for (int i = 0; i < capacity; i++)
		number[i] = 0xFFFFFFFF;
}

int BigDouble::highestNonZeroBit() {
	int i, j;
	for (i = capacity - 1; i > 0; i--)
		if (number[i] != 0)
			break;
	for (j = base - 1; j > 0; j--)
		if (((number[i] >> j) & 1) != 0)
			break;
	return (base * i + j);
}

int BigDouble::highestNonZeroBlock() {
	int i;
	for (i = capacity - 1; i > 0; i--)
		if (number[i] != 0)
			break;
	return i;
}