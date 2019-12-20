#ifndef SNGS_SNGS_H
#define SNGS_SNGS_H

#include "sbioinfo/sbam.h"
#include "sbioinfo/sbioseq.h"

namespace slib {
    namespace sbio {
        

        struct sngs_param {
			int32_t depth_bin, thread_count;
            
			bool parallel, ignore_dup;
            
			std::mutex mtx;
            
			sngs_param();
            ~sngs_param();
            
			sngs_param &operator=(const sngs_param &par);
            
            void set(const sobj &obj);
            sobj toObj();
        };
        
        
        class SNGSData {
        public:
            typedef Array<vararray> vararray2d;
            
        private:
            stpool _threads;
            std::mutex *_mtxs;
			//Array<SLock> _locks;
            
        public:
            int32_t ref_num, depth_bin;
            intarray ref_length, depth_size, uncovered;
            suinteger total_length, total_reads;
            bool target_seq;
            Array<sregion> target;
            integerarray read_count;
            double average_length, average_depth, covered_region;
            doublearray read_length;
            
			floatarray2d depth;

            vararray2d variants;
            intarray2d delidx, insidx, invidx;
            intarray3d trsidx, trinvidx;
            
        public:
            SNGSData();
            SNGSData(sngs_param *p);
            SNGSData(sngs_param *p, SBSeqList *list);
            SNGSData(sngs_param *p, SBamFile *bam);
            SNGSData(const char *path);
            ~SNGSData();
            
            void load(const char *path);
            void save(const char *path);
            
            void setNum(int32_t n);
            void setLength(int i, int32_t l);
            void setParam(sngs_param *p);
            
            void lock(int r, int v = 0);
            void unlock(int r, int v = 0);
            
            double depthAt(const sbpos &pos);
            /*
             double totalDepthIn(int32_t ref, int32_t pos, int32_t length);
             float depthIn(int32_t ref, int32_t pos, int32_t length);
             */
            void varindex(svariant_param *vp);
            void subtract(SNGSData &sum, svariant_param *vp);
            void tidy(svariant_param *vp);
            void integrate(SNGSData &sum, svariant_param *vp);
            
            void init();
        };
    }
}



#endif
