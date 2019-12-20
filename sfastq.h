#ifndef SNGS_SFASTQ_H
#define SNGS_SFASTQ_H

#include "sbioinfo/sbioutil.h"
#include "sbioinfo/sbioseq.h"

namespace slib {
	namespace sbio {
		class SFastq : public sbio::SBSeqList {
		public:
			SFastq();
			~SFastq();

			void importFq(const char* p);
			void exportFq(const char* p);
		};
	}
}

#endif
