#ifndef __SBIO_SVARYTIL_H__
#define __SBIO_SVARYTIL_H__

#include "sbioinfo/sbannot.h"

namespace slib {
	namespace sbio {

		class SVariant;
		class SVarList;

		extern String vtype(sushort i);
		extern sushort vtypeIdx(const char* s);
		extern String vsite(sushort i);
		extern subyte vsiteIdx(const char* s);
		extern String mtype(sushort i);
		extern subyte mtypeIdx(const char* s);
		

		struct smut_gene {
			String gene;
			sushort region;

			smut_gene();
			~smut_gene();

		};

		struct svar_annot {
			



		};



		class SVarIO {
		public:
			static const stringarray VCF_TABLE_COLUMNS;

			static void loadTxt(sio::SFile& file, SVarList* list, SBSeqList* ref);
			static void loadTCV(sio::SFile& file, SVarList* list, SBSeqList* ref);
			static void loadVCF(sio::SFile& file, SVarList* list, SBSeqList* ref, const sattribute& converter = {});
			static void loadJSON(sio::SFile& file, SVarList* list);
			static void BEDtoVar(String& str, SVariant* var);
			static void VCtoVar(String& str, SVariant* var);
			
			static void saveTxt(sio::SFile& file, SVarList* list);
			static void saveTCV(sio::SFile& file, SVarList* list, const stringarray& col = {""});
			static void saveVCF(sio::SFile& file, SVarList* list, SBSeqList* ref);
			static void saveJSON(sio::SFile& file, SVarList* list);
			static void toBED(String& str, SVariant* var);
			static void toVC(String& str, SVariant* var);

		};

		class SVarAnnotator {
		private:
			SBSeqList* _ref;
			SBAnnotDB *_db;

		public:
			SVarAnnotator(SBSeqList* ref = nullptr, SBAnnotDB* db = nullptr);
			~SVarAnnotator();

			void analyze(SVariant* var);
			void setReference(SBSeqList* ref);
			void setDB(SBAnnotDB* db);
		};

		class SVarFilter {
		private:


		public:
			SVarFilter();
			~SVarFilter();

			void checkRepeat(SVarList *list, SBSeqList* ref);
			void filter(SVarList* list);
			//void sortBy();
			bool isAvailable(SVariant* var);

		};


	}
}

#endif