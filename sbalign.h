#ifndef SBIO_SBALIGN_H
#define SBIO_SBALIGN_H

#include "sobj.h"
#include "sbioinfo/sbioseq.h"

namespace slib {
    namespace sbio {
        constexpr sint DEFAULT_ALIGN_LENGTH = 20;
        constexpr sbyte DEFAULT_PM_SCORE = 2;
        constexpr sbyte DEFAULT_AM_SCORE = 1;
        constexpr sbyte DEFAULT_MM_SCORE = -1;
        constexpr sbyte DEFAULT_GAP_SCORE = -2;
        constexpr sbyte DEFAULT_GAP2_SCORE = -1;
        
        struct salignment_param {
			sbyte pm_score, am_score, mm_score, gap_score, gap2_score;
            sushort seq_type;
            sshort align_length;
            String score_name;
            matb score_table;
            matb compare_table;
            
            salignment_param(uint8_t t = MISC_SEQ);
            ~salignment_param();
            
            salignment_param &operator=(const salignment_param &par);
            
            void set(const sobj &obj);
            sobj toObj();
            void makeTable();
        };
        
        struct scigar {
            static constexpr uint8_t MATCH = 0;
            static constexpr uint8_t INSERTION = 1;
            static constexpr uint8_t DELETION = 2;
            static constexpr uint8_t SKIP = 3;
            static constexpr uint8_t SCLIP = 4;
            static constexpr uint8_t HCLIP = 5;
            static constexpr uint8_t PADDING = 6;
            static constexpr uint8_t PMATCH = 7;
            static constexpr uint8_t MMATCH = 8;
            static const char *CIGAR_STRING;
            static const Map<char, uint8_t> CIGAR_INDEX;
            
            uint8_t option;
            int32_t length;
            
            scigar();
            scigar(uint8_t o, int32_t l);
            scigar(uint32_t i);
            scigar(const scigar &c);
            
            scigar & operator = (const scigar &c);
            bool operator < (const scigar &c) const;
            bool operator == (const scigar &c) const;
            bool operator != (const scigar &c) const;
        };
        
        class SCigarArray : public BiArray<scigar> {
        public:
            typedef BiArray<scigar> cigarray;
            
        public:
            SCigarArray();
            SCigarArray(size_t n);
            SCigarArray(const char *s);
            SCigarArray(const scigar &c);
            SCigarArray(int32_t n, uint32_t *c);
            SCigarArray(std::initializer_list<scigar> li);
            SCigarArray(const SCigarArray &c);
            ~SCigarArray();
            
            SCigarArray &operator = (const SCigarArray &c);
            
            void add(const scigar &c);
            void append(const SCigarArray &c);
            void put(const scigar &c);
            void pile(const SCigarArray &c);
            void set(int n, uint32_t *c);
            void reverse();
            
            size_t countRef(size_t beg = 0, size_t len = -1) const;
            size_t countQue(size_t beg = 0, size_t len = -1) const;
            size_t countCigar(int8_t op);
            
            std::string toString() const;
            
            bool operator==(const SCigarArray &array) const;
        };
        
        struct salign {
            sbpos ref;
            srange aligned;
            int32_t score;
            SCigarArray cigars;
            
            salign();
            salign(const sbpos &pos, const srange &range);
            salign(const salign &align);
            ~salign();
            
            salign &operator = (const salign &align);
            
            bool operator < (const salign &align) const;
            bool operator == (const salign &align) const;
            
            void scoring(salignment_param *par);
            void init();
            
            String alref(const String &ref);
            String match();
            String consensus(const String &ref, const String &que);
            String alque(const String &que);
        };
        
#define alignarray Array<salign>
        
        class SAlignment {
            
        private:
            salignment_param *_par;
            bytearray _path;
            intarray _score, _score2, _maxcol, _maxrow;
            int32_t _scr[4];
            
        public:
            intarray scores;
            SCigarArray cigars;
            
        public:
            SAlignment();
            SAlignment(salignment_param *p);
            ~SAlignment();
            
            void align(uint8_t *ref, size_t rlen, uint8_t *que, size_t qlen);
            void ralign(uint8_t *ref, size_t rlen, uint8_t *que, size_t qlen);
            void set(salignment_param *p);
            void init();
            void reset();
        };
        
        extern void calcScore(salign *al, matb *comparer, bytearray *score_list);
    }
}

#endif