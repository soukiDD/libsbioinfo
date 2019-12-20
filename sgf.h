#ifndef SBIO_SGF_H
#define SBIO_SGF_H

#include "sbioinfo/sbioutil.h"
#include "sbioinfo/sbannot.h"

namespace slib {
    using namespace sio;
    
    namespace sbio {
        
        struct sgff : public sbpos {
            String seqid, source, type, strand;
            sbyte phase;
            float score;
            sattribute attribute;
            
            sgff();
            sgff(const String &row);
            sgff(const sgff &g);
            ~sgff();
            
            void init();
            void set(const String &row);
        };
        
        class SGFFFile {
        private:
            SFile _file;
            String _row;
            bool _alias;
            
        protected:
            sgff gff_data;
            Array<sgff> gff_datas;
            
        public:
            SGFFFile();
            SGFFFile(const char *path, bool l);
            ~SGFFFile();
            
            void open(const char *path);
            void load(const char *path);
            
            const sgff &data() const;
            const Array<sgff> &list() const;
            
            bool next();
            bool loaded();
        };
    }
}

#endif
