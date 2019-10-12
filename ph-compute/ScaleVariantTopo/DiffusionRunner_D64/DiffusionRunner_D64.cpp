// -------------------------------------------------------------------
// $ProjectName     : DiffusionRunner_D64
// $FileName        : DiffusionRunner_D64.cpp$
// $Programmer      : Tran Quoc Hoan$
// -------------------------------------------------------------------

#include "stdafx.h"
#include "DiffusionRunner_D64.h"
#include "NpMatUtils.h"
#include "../TopoUtils_D64/StringUtils.h"
#include "../RipserLib_D64/RipserComputeUtils.h"

#include <boost/filesystem.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/foreach.hpp>
#include "boost/program_options.hpp"
#include <fstream>

#include <ppl.h>
#include <ppltasks.h>

using namespace NRipserComputeUtils;
using namespace NNpMatUtils;

std::vector<std::wstring> ReadMultiLinesFromFile(const std::wstring& filepath) {
	std::vector<std::wstring> rstr;
	std::ifstream read_op;
	read_op.open(filepath);
	if (!read_op.good()) {
		return rstr;
	}
	else {
		size_t index = 0;
		while (read_op.good()) {
			std::string line;
			std::getline(read_op, line);
			if (line.empty() == false)
				rstr.push_back(NStringUtil::_s2w(line));
		}
	}
	read_op.close();
	return rstr;
}

bool WriteDiffusionBarcodesToFile(std::vector<std::pair<size_t, RipsComputePrmPtr>> diffusion_barcodes,
	std::wstring out_dir, size_t maxdim, std::wstring basename) {
	namespace io = boost::iostreams;
	// Make output directory if not exist
	boost::filesystem::path outpath(out_dir);
	if (!boost::filesystem::is_directory(outpath)) {
		boost::filesystem::create_directory(outpath);
	}
	for (size_t dim = 0; dim <= maxdim; ++dim) {
		std::wstring outfile = out_dir + L"/diffusion_barcode_" + basename + L"_dim_" + std::to_wstring(dim) + L".txt";
		io::stream_buffer<io::file_sink> buf(NStringUtil::_w2s(outfile));
		std::ostream out(&buf);

		for (auto dbar : diffusion_barcodes) {
			size_t tau = dbar.first;
			RipsComputePrmPtr prm = dbar.second;
			if (!prm) continue;
			OutputPrmPtr output_prm = prm->output_prm;
			if (!output_prm) continue;
			auto bvec = output_prm->Barcodes();
			if (bvec.size() <= dim) continue;
			for (auto bar : bvec[dim]->barcodes()) {
				out << bar->birth << " " << bar->death << " " << tau << " " << std::endl;
			}
		}
        std::cout << "Saved diffusion barcodes to file: " << NStringUtil::_w2s(outfile) << std::endl;
	}
	return true;
}

int main(int argc, char** argv)
{
	using namespace boost::program_options;
	namespace fs = boost::filesystem;
	int nRetCode = 0;
	options_description description("Persistence runner for diffusion distance");
	description.add_options()
		("debug,x", value<bool>()->default_value(false), "Print Debug Information")
        ("norm", value<size_t>()->default_value(0), "Normalize the distance before calculating")
        ("nthreads", value<int>()->default_value(-1), "Number of threads to use, -1: use all as possible")
		("modulus,m", value<coefficient_t>()->default_value(2), "Compute homology with coefficients in the prime field Z/<p>Z")
		("maxdim,d", value<index_t>()->default_value(0), "Compute persistent homology up to dimension <k>")
		("thres", value<value_t>()->default_value(std::numeric_limits<value_t>::max()), "Compute Rips complexes up to diameter <t>")
		("outdir,o", value<std::string>()->default_value("output"), "Output directory")
		("input,i", value<std::string>()->default_value("network"), "Input folder")
		("taumax,t", value<size_t>()->default_value(0), "Compute with diffusion time up to taumax")
		("nbegin,b", value<size_t>()->default_value(0), "Begin index")
		("nend,e", value<size_t>()->default_value(0), "End index")
        ("skiptime,s", value<size_t>()->default_value(0), "Stop and continue for slow compute of Rips")
		("help,H", "Help: Usage DiffusionRunner_D64 [options]")
		("version,v", "v1.0")
		;
	variables_map vm;
	store(parse_command_line(argc, argv, description), vm);
	notify(vm);

	if (vm.count("help")) {
		std::cout << description << std::endl;
		return nRetCode;
	}
    size_t norm = vm["norm"].as<size_t>();
	size_t taumax = vm["taumax"].as<size_t>();
	size_t nbegin = vm["nbegin"].as<size_t>();
	size_t nend   = vm["nend"].as<size_t>();
    size_t skiptime = vm["skiptime"].as<size_t>();
    if (skiptime > 0) {
        std::cout << "Enable skiptime: " << skiptime << std::endl;
    }

    auto nthreads = vm["nthreads"].as<int>();
    if (nthreads <= 0) {
        nthreads = std::thread::hardware_concurrency();
    }
    std::cout << "Number of threads (= " << nthreads << ")" << std::endl;

    auto input_str = vm["input"].as<std::string>();
    std::cout << "Input: " << input_str.c_str() << std::endl;

	RipsComputePrmPtr prm(new RipsComputePrm());
    prm->skiptime = skiptime;
	RipsPrmPtr rip_prm = prm->rip_prm;
	InputPrmPtr input_prm = prm->input_prm;
	OutputPrmPtr output_prm = prm->output_prm;

	auto maxdim = vm["maxdim"].as<index_t>();
	auto out_dir = NStringUtil::_s2w(vm["outdir"].as<std::string>());
    std::cout << "Output: " << NStringUtil::_w2s(out_dir) << std::endl;

	rip_prm->modulus = vm["modulus"].as<coefficient_t>();
    if (false == NRipserUtils::isPrime(rip_prm->modulus)) {
        std::cout << "Modulus = (" << rip_prm->modulus << ") must be a prime number" << std::endl;
        return false;
    }
	rip_prm->dim_max = maxdim;
	rip_prm->diameter_max = vm["thres"].as<value_t>();

	output_prm->write_mode = INOUT_MODE::INNER_MODE;
	output_prm->out_dir = out_dir;
	input_prm->read_mode = INOUT_MODE::INNER_MODE;
	input_prm->format = LOWER_DISTANCE_MATRIX;
    
	// read all input distance matrices
	std::vector<fs::path> input_files;
	fs::path input_path(input_str);
	std::vector<std::wstring> rstr;

	if (fs::is_directory(input_path)) {
		// loop for all files in folder
		fs::directory_iterator it(input_path), eod;
		BOOST_FOREACH(fs::path const &p, std::make_pair(it, eod)) {
			if (fs::is_regular_file(p) && (fs::extension(p) == ".npy")) {
				input_files.push_back(p);
			}
		}
	}
	else if (fs::is_regular_file(input_path) && (fs::extension(input_path) == ".npy")) {
		input_files.push_back(input_path);
	}
	else {
        std::cout << "Error: no input files" << std::endl;
		return nRetCode;
	}

    size_t n_end_files = input_files.size();
    if (nend > 0) n_end_files = std::min(n_end_files, nend);
    std::cout << "Number of input files to process: " << n_end_files << std::endl;

	for (size_t i = nbegin; i < n_end_files; ++i) {
		auto p = input_files[i];
		NpyArray narr;
        const std::string pathname = NStringUtil::_w2s(p.c_str());
        //std::cout << "Process file: " << pathname << std::endl;
        if (false == NNpMatUtils::NpyLoad(pathname, narr)) {
            std::cout << "Error: loading npy file" << std::endl;
            return false;
        }
            
		const size_t num_taus = (size_t)narr.shape[0];
		const size_t Nr = (size_t)narr.shape[1];
		const size_t Nc = (size_t)narr.shape[2];
		const size_t N = Nr * Nc;
        // std::cout << "Size of loaded npy file: " << num_taus << "," << Nr << "," << Nc << std::endl;

		FType* dat = reinterpret_cast<FType*>(narr.data);
		std::vector<std::pair<size_t, RipsComputePrmPtr>> diffusion_barcodes;
		std::vector<RipsComputePrmPtr> prm_vec;

        size_t ntaus_to_process = num_taus;
        if (taumax > 0) ntaus_to_process = std::min(num_taus, taumax);

		for (size_t n = 0; n < ntaus_to_process; ++n) {
			RipsComputePrmPtr prm_f(new RipsComputePrm());
			prm_f->ParamCopy(prm);
            // Compute the maximum distance
            FType maxdist = dat[n*N];
            for (size_t i = 0; i < Nr; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    const FType tmp = dat[n*N + i*Nc + j];
                    if (tmp > maxdist) maxdist = tmp;
                }
            }
            if (maxdist <= 0) continue;
			// Make distance
			std::vector<value_t> distvec;
			for (size_t i = 0; i < Nr; ++i)
                for (size_t j = 0; j < i; ++j) {
                    FType tmp = dat[n*N + i*Nc + j];
                    if (norm > 0) tmp = tmp / maxdist;
                    distvec.push_back(tmp);
                }
			prm_f->input_prm->distances = distvec;
			prm_vec.push_back(prm_f);
			diffusion_barcodes.push_back(std::make_pair(n, prm_f));
		}
		if (!prm_vec.empty()) {
			ComputeRipPHMultiFiles(nthreads, prm_vec);
		}
		std::wstring basename = p.stem().c_str();
		WriteDiffusionBarcodesToFile(diffusion_barcodes, out_dir, maxdim, basename);
	}
	return nRetCode;
}