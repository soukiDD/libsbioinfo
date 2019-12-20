#ifndef SBIO_SBIOSEQ_H
#define SBIO_SBIOSEQ_H

#include "sio/sfile.h"
#include "sbioinfo/sbsequtil.h"
#include "sbioinfo/sbseqio.h"

namespace slib {
    namespace sbio {
        
        class SBioSeq;
        class SBSeqList;
        
        #define sbseq sptr<SBioSeq>
        
        constexpr size_t MAX_FIND_SEQ_LENGTH = 1<<14;
        
        class SBioSeq : public ubytearray {
			friend SBSeqIO;
            friend SBSeqList;
            
        protected:
            
        public:
            typedef Array<sbseq_annot> annot_vec;
            
        protected:
            sushort _type;
            sint _length;
            String _name;
            sregion _mask;
            annot_vec _annotation;
            SEQ_CONVERTER _dec, _enc;
            
        private:
            void _init();
            void setType();
            
        public:
            SBioSeq();
            SBioSeq(sushort t, const char *name = nullptr, const char *seq = nullptr);
            SBioSeq(const SBioSeq &seq);
            ~SBioSeq();
            
            SBioSeq &operator=(const SBioSeq &seq);
            
            void load(sushort t, const char *path);
            void save(const char *path);
            
            bool isAlias() const;
            bool isCircular() const;
            bool isMasked() const;
            bool isScored() const;
            bool isAnnotated() const;
            sushort type() const;
            subyte seqtype() const;
            subyte compress() const;
            
            const sint &length() const;
            const String &name() const;
            const sregion &mask() const;
            
            void setSeq(const char *seq);
            void setName(const char *name);
            void setLength(const size_t &l, bool alias = false);
            void addMask(const srange &range);
            void removeMask(const srange &range);
            void clearAll();
            
            void encode(const char *seq, size_t off = 0, size_t len = -1, bool dir = false);
            void decode(char *seq, size_t off = 0, size_t len = -1, bool dir = false);
            void recode(subyte c, ubytearray &seq, size_t off = 0, size_t len = -1, bool dir = false);
            
            void convert(sushort t);
            SBioSeq subseq(const sbpos &pos);
            SBioSeq subseq(size_t off, size_t len, bool dir = false);
            String raw(const sbpos &pos, bool mask = false) const;
            String raw(size_t off = 0, size_t len = -1, bool dir = false, bool mask = false) const;
            
            void append(const char *seq);
            void append(const SBioSeq &seq);
            void insert(size_t idx, const char *seq);
            void insert(size_t idx, const SBioSeq &seq);
            
            void circulate();
            void complement();
            void splice(const sregion &region);
            void transcribe();
            void rtranscribe();
            void translate(const CODON_TABLE &code = DEFAULT_CODON);
            
            SBioSeq transcript();
            SBioSeq rtranscript();
            SBioSeq translated(const CODON_TABLE &code = DEFAULT_CODON);
        };
        
        #define sdna(X) sbseq(COMPRESS1|DNA_SEQ, nullptr, (X))
        #define srna(X) sbseq(COMPRESS1|RNA_SEQ, nullptr, (X))
        #define saa(X) sbseq(COMPRESS1|AA_SEQ, nullptr, (X))
        
        class SBSeqList : public slib::Array<sbseq> {
            typedef uintegerarray offset_vec;
            
        protected:
            sindex _index;
            offset_vec _offset;
            sio::SFile _file;
            
        public:
            SBSeqList();
            SBSeqList(const char *s, bool l = true);
            virtual ~SBSeqList();
            
            const sindex &nameIdx() const;
            size_t seqIdx(const char *name);
            suinteger total() const;
            
            void load(const char *path);
            void save(const char *path);
            void index(const char *path);
            
            void readSeq(int idx);
            void writeSeq(int idx);
            
            SBioSeq getSeq(sushort type, const sbpos &pos);
            SBioSeq getSeq(sushort type, int idx, size_t pos, size_t len, bool dir = false);
            
            sio::SFile &file();
        };
        
        class SFasta : public SBSeqList {
        private:
            size_t row_count;
        public:
            SFasta();
            ~SFasta();
            
            void importFa(const char *path, sushort type = 0, size_t init = 256);
            void exportFa(const char *path);
            void indexFa(const char *path, sushort type = 0);
            void importSeq(int idx, sushort type = 0);
            void exportSeq(int idx);
        };
    }
}

#endif
