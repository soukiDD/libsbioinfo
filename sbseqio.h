#pragma once

#include "sobj.h"

namespace slib {
	namespace sbio {

		constexpr size_t FASTA_ROW_CHAR = 60;

		class SBioSeq;
		class SBSeqList;

		struct abidir {
			sint number;
			sshort element_type;
			sshort element_size;
			sint num_elements;
			sint data_size;
			sint data_offset;
			sint data_handle;

			abidir();
			abidir(const abidir& dir);
			~abidir();

			abidir& operator=(const abidir& dir);
			bool operator<(const abidir& dir) const;
			bool operator==(const abidir& dir) const;
		};

		class SBSeqIO {

		public:
			SBSeqIO();
			~SBSeqIO();

			static void loadTXT(sio::SFile& file, SBioSeq* seq);
			static void loadABI(sio::SFile& file, SBioSeq* seq);
			static void loadGBK(sio::SFile& file, SBioSeq* seq);
			static void loadFASTA(sushort type, sio::SFile& file, SBioSeq* seq);
			static void loadFASTA(sushort type, sio::SFile& file, SBSeqList* list);
			static void saveTXT(sio::SFile& file, SBioSeq* seq);
			static void saveGBK(sio::SFile& file, SBioSeq* seq);
			static void saveFASTA(sio::SFile& file, SBioSeq* seq);
			static void saveFASTA(sio::SFile& file, SBSeqList* list);
		};

	}
}