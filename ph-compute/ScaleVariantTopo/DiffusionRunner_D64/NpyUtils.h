#pragma once
#include "stdafx.h"
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <cstdio>
#include <typeinfo>
#include <iostream>
#include <cassert>
#include <map>
#include <boost/algorithm/string/join.hpp>
#include <boost/range/adaptor/transformed.hpp>

namespace NNpyUtils {
	struct NpyArray {
		char* data;
		std::vector<size_t> shape;
		size_t word_size;
		bool fortran_order;
		const size_t count() const {
			size_t count = 1;
			for (size_t s : shape) count *= s;
			return count;
		}
		std::string shape_to_str() {
			using boost::algorithm::join;
			using boost::adaptors::transformed;

			auto tostr = static_cast<std::string(*)(size_t)>(std::to_string);
			return join(shape | transformed(tostr), ",");
		}
		void destruct() { delete[] data; }
	};

	struct npz_t : public std::map<std::string, NpyArray> {
		void destruct() {
			npz_t::iterator it = this->begin();
			for (; it != this->end(); ++it) (*it).second.destruct();
		}
		bool key_check(std::string key) {
			return (this->count(key) > 0);
		}
	};

	template<typename T>
	std::string ToString(T i, int pad = 0, char padval = ' ') {
		std::stringstream s;
		s << i;
		return s.str();
	}

	inline char BigEndianTest() {
		unsigned char x[] = { 1, 0 };
		short y = *(short*)x;
		return y == 1 ? '<' : '>';
	}

	inline char MapType(const std::type_info& t) {
		if (t == typeid(float)) return 'f';
		if (t == typeid(double)) return 'f';
		if (t == typeid(long double)) return 'f';

		if (t == typeid(int)) return 'i';
		if (t == typeid(char)) return 'i';
		if (t == typeid(short)) return 'i';
		if (t == typeid(long)) return 'i';
		if (t == typeid(long long)) return 'i';

		if (t == typeid(unsigned char)) return 'u';
		if (t == typeid(unsigned short)) return 'u';
		if (t == typeid(unsigned long)) return 'u';
		if (t == typeid(unsigned long long)) return 'u';
		if (t == typeid(unsigned int)) return 'u';

		if (t == typeid(bool)) return 'b';

		else return '?';
	}

	template<typename T>
	std::vector<char> CreateNpyHeader(const T* data, const std::vector<size_t> shape, const size_t ndims);

	bool ParseNpyHeader(FILE* fp, size_t& word_size, std::vector<size_t>& shape, size_t& ndims, bool& fortran_order);
	bool LoadNpyFile(FILE* fp, NpyArray& arr);
	bool NpzLoad(std::string fname, npz_t& arrays);
	bool NpzLoad(std::string fname, std::string varname, NpyArray& array);
    bool NpyLoad(std::string fname, NpyArray& arr);

	template<typename T>
	std::vector<char>& operator+=(std::vector<char>& lhs, const T rhs);

	template<>
	std::vector<char>& operator+=(std::vector<char>& lhs, const std::string rhs);

	template<>
	std::vector<char>& operator+=(std::vector<char>& lhs, const char* rhs);

	template<typename T>
	bool NpySave(std::string fname, const T* data, const std::vector<size_t> shape, const size_t ndims, std::string mode = "w");

	bool NpySaveF(std::string fname, const float* data, const std::vector<size_t> shape, const size_t ndims, std::string mode = "w");
}

