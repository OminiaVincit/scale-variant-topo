#include "stdafx.h"
#include "NpyUtils.h"

namespace NNpyUtils {
	template<typename T>
	std::vector<char> NNpyUtils::CreateNpyHeader(const T* data, const std::vector<size_t> shape, const size_t ndims) {
		std::vector<char> dict;
		dict += "{'descr': '";
		dict += BigEndianTest();
		dict += MapType(typeid(T));
		dict += ToString(sizeof(T));
		dict += "', 'fortran_order': False, 'shape': (";
		dict += ToString(shape[0]);
		for (size_t i = 1; i < ndims; i++) {
			dict += ", ";
			dict += ToString(shape[i]);
		}
		if (ndims == 1) dict += ",";
		dict += "), }";
		//pad with spaces so that preamble+dict is modulo 16 bytes. preamble is 10 bytes. dict needs to end with \n
		int remainder = 16 - (10 + dict.size()) % 16;
		dict.insert(dict.end(), remainder, ' ');
		dict.back() = '\n';

		std::vector<char> header;
		header += (char)0x93;
		header += "NUMPY";
		header += (char)0x01; //major version of numpy format
		header += (char)0x00; //minor version of numpy format
		header += (unsigned short)dict.size();
		header.insert(header.end(), dict.begin(), dict.end());

		return header;
	}

	bool NNpyUtils::ParseNpyHeader(FILE* fp, size_t& word_size, std::vector<size_t>& shape, size_t& ndims, bool& fortran_order) {
		char buffer[256];
		size_t res = fread(buffer, sizeof(char), 11, fp);
        
		if (res != 11)
			throw std::runtime_error("ParseNpyHeader: failed fread");
		std::string header = fgets(buffer, 256, fp);
        if (header.empty() || header[header.size() - 1] != '\n') {
            std::cout << "Error: header in ParseNpyHeader: " << header << std::endl;
            return false;
        }

		size_t loc1, loc2;

		//fortran order
		loc1 = header.find("fortran_order") + 16;
		fortran_order = (header.substr(loc1, 5) == "True" ? true : false);

		//shape
		{
			loc1 = header.find("(");
			loc2 = header.find(")");
			std::string str_shape = header.substr(loc1 + 1, loc2 - loc1 - 1);
			if (str_shape[str_shape.size() - 1] == ',') 
                ndims = 1;
			else
                ndims = std::count(str_shape.begin(), str_shape.end(), ',') + 1;
			std::vector<size_t> newshape(ndims);
			for (size_t i = 0; i < ndims; ++i) {
				loc1 = str_shape.find(",");
				newshape[i] = atoi(str_shape.substr(0, loc1).c_str());
				str_shape = str_shape.substr(loc1 + 1);
			}
			shape = newshape;
		}

		//endian, word size, data type
		//byte order code | stands for not applicable. 
		//not sure when this applies except for byte array
		loc1 = header.find("descr") + 9;
		bool littleEndian = (header[loc1] == '<' || header[loc1] == '|' ? true : false);
		if(false == littleEndian)
            return false;

		//char type = header[loc1+1];
		//assert(type == map_type(T));

		std::string str_ws = header.substr(loc1 + 2);
		loc2 = str_ws.find("'");
		word_size = atoi(str_ws.substr(0, loc2).c_str());
		return true;
	}

	bool NNpyUtils::LoadNpyFile(FILE* fp, NpyArray& arr) {
		std::vector<size_t> shape;
		size_t ndims, word_size;
		bool fortran_order;
		if (false == ParseNpyHeader(fp, word_size, shape, ndims, fortran_order))
            return false;
		size_t size = 1;
		for (size_t i = 0; i < ndims; ++i) size *= shape[i];

		arr.word_size = word_size;
		arr.shape = shape;
		arr.data = new char[size*word_size];
		arr.fortran_order = fortran_order;
		size_t nread = fread(arr.data, word_size, size, fp);
		if (nread != size) {
			return false;
		}
		return true;
	}

	bool NNpyUtils::NpzLoad(std::string fname, npz_t& arrays) {
		FILE* fp;
		fopen_s(&fp, fname.c_str(), "rb");

        if (!fp) return false;

		while (1) {
			std::vector<char> local_header(30);
			size_t headerres = fread(&local_header[0], sizeof(char), 30, fp);
			if(headerres != 30)
                return false;

			//if we've reached the global header, stop reading
			if (local_header[2] != 0x03 || local_header[3] != 0x04) break;

			//read in the variable name
			unsigned short name_len = *(unsigned short*)&local_header[26];
			std::string varname(name_len, ' ');
			size_t vname_res = fread(&varname[0], sizeof(char), name_len, fp);

			if (vname_res != name_len)
                return false;

			//erase the lagging .npy        
			varname.erase(varname.end() - 4, varname.end());

			//read in the extra field
			unsigned short extra_field_len = *(unsigned short*)&local_header[28];
			if (extra_field_len > 0) {
				std::vector<char> buff(extra_field_len);
				size_t efield_res = fread(&buff[0], sizeof(char), extra_field_len, fp);
				if (efield_res != extra_field_len)
                    return false;
			}
			if(false == LoadNpyFile(fp, arrays[varname]))
                return false;
		}

		fclose(fp);
		return true;
	}

	bool NNpyUtils::NpzLoad(std::string fname, std::string varname, NpyArray& array) {
		FILE* fp;
		fopen_s(&fp, fname.c_str(), "rb");

		if(!fp) return false;

		while (1) {
			std::vector<char> local_header(30);
			size_t header_res = fread(&local_header[0], sizeof(char), 30, fp);
            if (header_res != 30)
                return false;

			//if we've reached the global header, stop reading
			if (local_header[2] != 0x03 || local_header[3] != 0x04) break;

			//read in the variable name
			unsigned short name_len = *(unsigned short*)&local_header[26];
			std::string vname(name_len, ' ');
			size_t vname_res = fread(&vname[0], sizeof(char), name_len, fp);

            if (vname_res != name_len)
                return false;
			vname.erase(vname.end() - 4, vname.end()); 
			unsigned short extra_field_len = *(unsigned short*)&local_header[28];
			fseek(fp, extra_field_len, SEEK_CUR); //skip past the extra field

			if (vname == varname) {
				LoadNpyFile(fp, array);
				fclose(fp);
				return true;
			}
			else {
				//skip past the data
				size_t size = *(size_t*)&local_header[22];
				fseek(fp, static_cast<long>(size), SEEK_CUR);
			}
		}

		fclose(fp);
		return false;
	}

    bool NNpyUtils::NpyLoad(std::string fname, NpyArray& arr) {
        
		FILE* fp;
		fopen_s(&fp, fname.c_str(), "rb");

        if (!fp) {
            return false;
        }
		LoadNpyFile(fp, arr);
		fclose(fp);
		return true;
	}

	template<typename T>
	bool NNpyUtils::NpySave(std::string fname, const T* data, const std::vector<size_t> shape, const size_t ndims, std::string mode) {
		FILE* fp = NULL;
		if (mode == "a") fopen_s(&fp, fname.c_str(), "r+b");
		if (fp) {
			//file exists. we need to append to it. read the header, modify the array size
			size_t word_size, tmp_dims;
			//size_t* tmp_shape = 0;
			std::vector<size_t> tmp_shape;
			bool fortran_order;
			ParseNpyHeader(fp, word_size, tmp_shape, tmp_dims, fortran_order);
			if(fortran_order)
                return false;

			if (word_size != sizeof(T)) {
                return false;
			}
			if (tmp_dims != ndims) {
                return false;
			}

			for (size_t i = 1; i < ndims; i++) {
				if (shape[i] != tmp_shape[i]) {
                    return false;
				}
			}
			tmp_shape[0] += shape[0];

			fseek(fp, 0, SEEK_SET);
			std::vector<char> header = CreateNpyHeader(data, tmp_shape, ndims);
			fwrite(&header[0], sizeof(char), header.size(), fp);
			fseek(fp, 0, SEEK_END);

			//delete[] tmp_shape;
		}
		else {
			fopen_s(&fp, fname.c_str(), "wb");
			std::vector<char> header = CreateNpyHeader(data, shape, ndims);
			fwrite(&header[0], sizeof(char), header.size(), fp);
		}

		size_t nels = 1;
		for (size_t i = 0; i < ndims; i++) nels *= shape[i];

		fwrite(data, sizeof(T), nels, fp);
		fclose(fp);
		return true;
	}

	bool NNpyUtils::NpySaveF(std::string fname, const float* data, const std::vector<size_t> shape, const size_t ndims, std::string mode) {
		return NpySave(fname, data, shape, ndims, mode);
	}

	template<typename T>
	std::vector<char>& NNpyUtils::operator+=(std::vector<char>& lhs, const T rhs) {
		//write in little endian
		for (char byte = 0; byte < sizeof(T); byte++) {
			char val = *((char*)&rhs + byte);
			lhs.push_back(val);
		}
		return lhs;
	}

	template<>
	std::vector<char>& NNpyUtils::operator+=(std::vector<char>& lhs, const std::string rhs) {
		lhs.insert(lhs.end(), rhs.begin(), rhs.end());
		return lhs;
	}

	template<>
	std::vector<char>& NNpyUtils::operator+=(std::vector<char>& lhs, const char* rhs) {
		//write in little endian
		size_t len = strlen(rhs);
		lhs.reserve(len);
		for (size_t byte = 0; byte < len; byte++) {
			lhs.push_back(rhs[byte]);
		}
		return lhs;
	}
}