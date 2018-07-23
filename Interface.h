#pragma once

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<ctime>

using namespace std;

const int BASE = 32;
const int ARRAY_LENGTH = 64;
const int DOUBLE_ARRAY_LENGTH = 2 * ARRAY_LENGTH;

const int FIELD_EXTENSION = 293;
const int VECTOR_LENGTH = 10; 
const int DENSITY = 5;

const int MILLER_RABIN_PRIMALITY_TEST_CONST = 15;

class BigInteger;
class BigDouble;
class PolynomialPB;
class PolynomialNB;

class BigInteger {
private:
	static const int base = BASE;
	static const int capacity = ARRAY_LENGTH;

	unsigned int* number;
public:
	friend class BigDouble;

	BigInteger();
	BigInteger(string n);
	BigInteger(unsigned int n);

	~BigInteger();

	string toString();

	void set(string n);
	void set(unsigned int n);
	void setBit(int n);

	void nullify();
	void maxNumber();

	int highestNonZeroBit();
	int highestNonZeroBlock();

	friend void copy(const BigInteger& scr, BigInteger& dst);
	friend void copy(const BigDouble& scr, BigInteger& dst);
	friend void copy(const BigInteger& scr, BigDouble& dst);

	friend int highestNonZeroBit(const BigInteger& x);
	friend int highestNonZeroBlock(const BigInteger& x);

	friend void add(const BigInteger& x, const BigInteger& y, BigInteger& r, unsigned int& carry);
	friend void add(const BigInteger& x, const BigInteger& y, BigInteger& r);
	friend void add(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend void sub(const BigInteger& x, const BigInteger& y, BigInteger& r, unsigned int& borrow);
	friend void sub(const BigInteger& x, const BigInteger& y, BigInteger& r);

	friend void increment(const BigInteger& x, BigInteger& r);

	friend void mult(const BigInteger& x, const unsigned int y, BigInteger& r);
	friend void mult(const BigInteger& x, const BigInteger& y, BigInteger& r);
	friend void mult(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend void div(const BigInteger& x, const BigInteger& y, BigInteger& q, BigInteger& r);

	friend void mod(const BigInteger& x, const BigInteger& n, BigInteger& r);

	friend void power(const BigInteger& x, const BigInteger& e, BigInteger& r);
	friend void powerWindow(const BigInteger& x, const BigInteger& e, BigInteger& r);

	friend void pow(const PolynomialPB& x, const BigInteger& e, PolynomialPB& r);
	friend void pow(const PolynomialNB& x, const BigInteger& e, PolynomialNB& r);

	friend int cmp(const BigInteger& x, const BigInteger& y);

	friend void leftShift(const BigInteger& x, int n, BigInteger& r);
	friend void leftBlocksShift(const BigInteger& x, int n, BigInteger& r);
	friend void rightShift(const BigInteger& x, int n, BigInteger& r);
	friend void rightBlocksShift(const BigInteger& x, int n, BigInteger& r);

	friend void gcd(const BigInteger& x, const BigInteger& y, BigInteger& r);
	friend void lcp(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend void inverseElement(const BigInteger& x, const BigInteger& n, BigInteger& r);

	friend void modAdd(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r);
	friend void modSub(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r);

	friend void modMult(const BigInteger& x, const BigInteger& y, const BigInteger& n, BigInteger& r);
	friend void modSqr(const BigInteger& x, const BigInteger& n, BigInteger& r);

	friend void barrettReduction(const BigDouble& x, const BigInteger& n, const BigDouble& m, BigInteger& r);
	friend void barrettConst(const BigInteger& n, BigDouble& m);
	friend void modPowBarrett(const BigInteger& x, const BigInteger& e, const BigInteger& n, BigInteger& r);

	friend void binaryAlgorithm(const BigInteger& x, const BigInteger& y, BigInteger& r);

	friend void linearCongruencesSystem(const BigInteger* b, const BigInteger* n, int size, BigInteger& r);

	friend void random(const int blocks, BigInteger& r);
	friend void MillerRabinPrimalityTest(const BigInteger& x, int& r);
};

class BigDouble {
private:
	static const int base = BASE;
	static const int capacity = DOUBLE_ARRAY_LENGTH;

	unsigned int* number;
public:
	friend class BigInteger;

	BigDouble();
	BigDouble(unsigned int n);
	BigDouble(string n);

	~BigDouble();

	string toString();

	void set(unsigned int n);
	void set(string n);
	void setBit(int n);

	void nullify();
	void maxNumber();

	int highestNonZeroBit();
	int highestNonZeroBlock();

	friend void copy(const BigDouble& scr, BigDouble& dst);
	friend void copy(const BigDouble& scr, BigInteger& dst);
	friend void copy(const BigInteger& scr, BigDouble& dst);

	friend int highestNonZeroBit(const BigDouble& x);
	friend int highestNonZeroBlock(const BigDouble& x);

	friend void add(const BigDouble& x, const BigDouble& y, BigDouble& r);
	friend void add(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend void sub(const BigDouble& x, const BigDouble& y, BigDouble& r);

	friend void increment(const BigDouble& x, BigDouble& r);

	friend void mult(const BigDouble& x, const unsigned int y, BigDouble& r);
	friend void mult(const BigDouble& x, const BigDouble& y, BigDouble& r);
	friend void mult(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend void div(const BigDouble& x, const BigDouble& y, BigDouble& q, BigDouble& r);

	friend void mod(BigDouble& x, BigDouble& n, BigDouble& r);

	friend void lcp(const BigInteger& x, const BigInteger& y, BigDouble& r);

	friend int cmp(const BigDouble &x, const BigDouble& y);

	friend void leftShift(const BigDouble& x, int n, BigDouble& r);
	friend void leftBlocksShift(const BigDouble& x, int n, BigDouble& r);
	friend void rightShift(const BigDouble& x, int n, BigDouble& r);
	friend void rightBlocksShift(const BigDouble& x, int n, BigDouble& r);

	friend void barrettReduction(const BigDouble& x, const BigInteger& n, const BigDouble& m, BigInteger& r);
	friend void barrettConst(const BigInteger& n, BigDouble& m);
};

class PolynomialPB {
private:
	static const int fieldExtension = FIELD_EXTENSION;

	static int irreduciblePolynomial[DENSITY];

	static const int base = BASE;
	static const int length = VECTOR_LENGTH;

	unsigned int* vector;
public:
	PolynomialPB();
	PolynomialPB(string n);

	~PolynomialPB();

	void addIdentity();
	void multIdentity();

	void set_bin(string n);
	void set_hex(string n);

	string toString();
	string toStringHex();

	int degree();

	friend void copy(const PolynomialPB& scr, PolynomialPB& dst);

	friend void add(const PolynomialPB& x, const PolynomialPB& y, PolynomialPB& r);

	friend void mult(const PolynomialPB& x, const PolynomialPB& y, PolynomialPB& r);

	friend int trace(const PolynomialPB& x);

	friend void sqr(const PolynomialPB& x, PolynomialPB& r);
	friend void pow(const PolynomialPB& x, unsigned int e, PolynomialPB& r);
	friend void pow(const PolynomialPB& x, const BigInteger& e, PolynomialPB& r);

	friend void inverse(const PolynomialPB& x, PolynomialPB& r);
};

class PolynomialNB {
private:
	static const int fieldExtension = FIELD_EXTENSION;

	static bool status;
	static unsigned int multiplicativeMatrix[FIELD_EXTENSION][2][2];

	static const int base = BASE;
	static const int length = VECTOR_LENGTH;

	unsigned int* vector;

	void calculatingMultiplicativeMatrix();

public:
	PolynomialNB();
	PolynomialNB(string n);

	~PolynomialNB();

	void addIdentity();
	void multIdentity();

	void set_bin(string n);
	void set_hex(string n);

	string toString();
	string toStringHex();

	friend void copy(const PolynomialNB& scr, PolynomialNB& dst);

	friend void add(const PolynomialNB& x, const PolynomialNB& y, PolynomialNB& r);

	friend void mult(const PolynomialNB& x, const PolynomialNB& y, PolynomialNB& r);
	friend void tensorMult(const PolynomialNB& x, const PolynomialNB& y, unsigned int& r);

	friend int trace(PolynomialNB& x);

	friend void sqr(const PolynomialNB& x, PolynomialNB& r);
	friend void pow(const PolynomialNB& x, unsigned int e, PolynomialNB& r);
	friend void pow(const PolynomialNB& x,  const BigInteger& e, PolynomialNB& r);

	friend void inverse(const PolynomialNB& x, PolynomialNB& r);

	friend void leftRotation(const PolynomialNB& x, int n, PolynomialNB& r);
	friend void rightRotation(const PolynomialNB& x, int n, PolynomialNB& r);
};