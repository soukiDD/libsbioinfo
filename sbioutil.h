#ifndef SBIO_SBIOINFO_UTIL_H
#define SBIO_SBIOINFO_UTIL_H

#include "sbio/sbio.h"

namespace slib {
    namespace sbio {

		constexpr sushort UNAVAILABLE_FLAG = 0x7000;

        //Annotation category
        constexpr uint16_t MISC_FEATURE = 0x0000;
        constexpr uint16_t CHROMOSOMIC_FEATURE = 0x0001;
        constexpr uint16_t GENOMIC_FEATURE = 0x0002;
        constexpr uint16_t GENE_FEATURE = 0x0003;
        constexpr uint16_t TRANSPOSON_FEATURE = 0x0004;
        constexpr uint16_t NUCLEOTIDE_FEATURE = 0x0005;
        constexpr uint16_t TRANSCRIPT_FEATURE = 0x0006;
        constexpr uint16_t PROTEIN_FEATURE = 0x0007;
        constexpr uint16_t VARIATION_FEATURE = 0x0008;
        
        //Chromosome type
        constexpr uint16_t NUCLEAR_CHR = 0x1000;
        constexpr uint16_t AUTOSOME = 0x1001;
        constexpr uint16_t SEX_CHR = 0x1F00;
        
        constexpr uint16_t EXT_CHR = 0x2000;
        constexpr uint16_t MT_GENOME = 0x2001;
        constexpr uint16_t PLASMID = 0x2002;
        constexpr uint16_t PLASTID = 0x2003;
        
        constexpr uint16_t ART_CHR = 0x4000;
        constexpr uint16_t BAC = 0x4001;
        constexpr uint16_t YAC = 0x4002;
        
        constexpr uint16_t BALANCER_CHR = 0x4100;
        
        //Gene type
        constexpr uint16_t PROTEN_CODING = 0x1000;
        constexpr uint16_t NON_CODING = 0x2000;
        constexpr uint16_t PSEUDOGENE = 0x4000;
        constexpr uint16_t TRANSPOSON = 0x8000;
        
        //Transcript type
        constexpr uint16_t M_RNA = 0x1000;
        
        constexpr uint16_t NC_RNA = 0x2000;
        constexpr uint16_t T_RNA = 0x2001;
        constexpr uint16_t R_RNA = 0x2002;
        constexpr uint16_t AS_RNA = 0x2003;
        constexpr uint16_t LNC_RNA = 0x2004;
        constexpr uint16_t LINC_RNA = 0x2005;
        
        constexpr uint16_t SMALL_RNA = 0x2100;
        constexpr uint16_t MI_RNA = 0x2101;
        constexpr uint16_t SI_RNA = 0x2102;
        constexpr uint16_t SN_RNA = 0x2103;
        constexpr uint16_t SNO_RNA = 0x2104;
        constexpr uint16_t PI_RNA = 0x2105;
        constexpr uint16_t SC_RNA = 0x2106;
        constexpr uint16_t SCA_RNA = 0x2107;
        constexpr uint16_t Y_RNA = 0x2108;
        constexpr uint16_t VT_RNA = 0x2109;
        constexpr uint16_t P_RNA = 0x210A;
        constexpr uint16_t SL_RNA = 0x210B;
        
        constexpr uint16_t EX_RNA = 0x21A0;
        
        constexpr uint16_t RIBOZYME = 0x2200;
        constexpr uint16_t RNASE_P = 0x2201;
        constexpr uint16_t RN_MRP = 0x2202;
        
        constexpr uint16_t TERC = 0x22A0;
        
        constexpr uint16_t VIRAL_RNA = 0x2400;
        constexpr uint16_t VA_RNA = 0x2401;
        
        constexpr uint16_t SH_RNA = 0x2800;
        
        constexpr uint16_t CRISPR_RNA = 0x2810;
        constexpr uint16_t GUIDE_RNA = 0x2811;
        constexpr uint16_t TRACER_RNA = 0x2812;
        
        //Genetic region type
        constexpr uint16_t CDS = 0x0001;
        constexpr uint16_t UTR = 0x0002;
        constexpr uint16_t UTR5 = 0x0012;
        constexpr uint16_t UTR3 = 0x0022;
        constexpr uint16_t EXON = 0x0004;
        constexpr uint16_t INTRON = 0x0008;
        constexpr uint16_t SPLICE_SITE = 0x0048;
        constexpr uint16_t PROCESSED = 0x0080;
        
        //Genomic feature type
        constexpr uint16_t CENTROMERE = 0x0001;
        constexpr uint16_t TEROMERE = 0x0002;
        
        constexpr uint16_t SATELLITE = 0x0010;
        
        constexpr uint16_t BALANCED_SITE = 0x0080;
        
        constexpr uint16_t REGULATOR = 0x0101;
        constexpr uint16_t PROMOTER = 0x0102;
        constexpr uint16_t ENHANCER = 0x0103;
        constexpr uint16_t REPRESSOR = 0x0104;
        constexpr uint16_t INSULATOR = 0x0105;
        constexpr uint16_t TERMINATOR = 0x0106;
        constexpr uint16_t SILENCER = 0x0107;
        
        constexpr uint16_t OPERON = 0x010A;
        constexpr uint16_t OPERATOR = 0x010B;
        
        constexpr uint16_t TSS_SITE = 0x0111;
        constexpr uint16_t TST_SITE = 0x0112;
        constexpr uint16_t TATA_BOX = 0x0113;
        
        constexpr uint16_t SPLICE_LEADER = 0x01A0;
        
        constexpr uint16_t LOW_COMPLEX = 0x0200;
        constexpr uint16_t POLY_A_SITE = 0x0201;
        constexpr uint16_t AT_RICH = 0x0202;
        constexpr uint16_t GC_RICH = 0x0203;
        
        constexpr uint16_t REPEAT_SITE = 0x0210;
        constexpr uint16_t TANDEM_REPEAT = 0x0211;
        constexpr uint16_t PALINDROME = 0x0212;
        
        constexpr uint16_t BINDING_SITE = 0x0400;
        constexpr uint16_t TF_SITE = 0x0401;
        
        
        constexpr uint16_t HETERO_CHR = 0x0801;
        constexpr uint16_t IMPRINT_SITE = 0x0802;
        
        constexpr uint16_t DNA_METHYL_SITE = 0x0810;
        
        constexpr uint16_t HISTONE_SITE = 0x0820;
        constexpr uint16_t HIS_METHYL_SITE = 0x0821;
        constexpr uint16_t HIS_ACETHYL_SITE = 0x0822;
        constexpr uint16_t HIS_PHOS_SITE = 0x0823;
        constexpr uint16_t HIS_UBIQ_SITE = 0x0824;
        
        constexpr uint16_t RESTRICTION_SITE = 0x1000;
        
        constexpr uint16_t PAM_SITE = 0x1F00;
        constexpr uint16_t CAS9_PAM = 0x1F01;
        constexpr uint16_t CAS12_PAM = 0x1F02;
        constexpr uint16_t CAS13_PAM = 0x1F03;
        
        //Variant type
        constexpr uint8_t SNV = 0x01;
        constexpr uint8_t MNV = 0x02;
        constexpr uint8_t INSERTION = 0x04;
        constexpr uint8_t DELETION = 0x08;
        constexpr uint8_t DUPLICATION = 0x10;
        constexpr uint8_t MULTIPLICATION = 0x20;
        constexpr uint8_t INVERSION = 0x40;
        constexpr uint8_t TRANSLOCATION = 0x80;
        
        //Genetic mutation type
        constexpr uint16_t SILENT = 0x0000;
        constexpr uint16_t SUBSTITUTION = 0x1000;
        constexpr uint16_t MISSENSE = 0x1100;
        constexpr uint16_t NONSENSE = 0x1200;
        constexpr uint16_t INDEL = 0x2000;
        constexpr uint16_t IN_FRAME = 0x2400;
        constexpr uint16_t FRAME_SHIFT = 0x2800;
        constexpr uint16_t STRUCTURE_CHANGE = 0x4000;
        constexpr uint16_t GENE_FUSION = 0x8000;
        /*
         constexpr uint16_t DEL_MUT = 0x0010;
         constexpr uint16_t MULTI_MUT = 0x0020;
         constexpr uint16_t PARTIAL_MUT = 0x1000;
         constexpr uint16_t COMPLETE_MUT = 0x1000;
         
         constexpr uint16_t REARRANGE = 0x0100;
         constexpr uint16_t TANDEM_ARRANGE = 0x0400;
         constexpr uint16_t WRAP_ARRANGE = 0x0800;
         
         
         constexpr uint16_t IPSI_DIR = 0x1000;
         constexpr uint16_t CONTRA_DIR = 0x2000;
         */
        extern float scoreVal(float v);
        extern float phredVal(float v);
        
        extern size_t countBin(sorder &order, srange range);
        extern size_t getBin(srange range);
        extern void getBins(sizearray &bins, srange range);
        extern void getBins(sizearray &bins, const sregion &region);
         
        struct sbpos : public srange {
            int32_t idx;
            bool dir;
            
            sbpos();
            sbpos(int i, int b, int e, bool d = false);
            sbpos(const char *s, sindex *namei = nullptr);
            sbpos(const sbpos &p);
            sbpos &operator = (const sbpos &p);
            ~sbpos();
            
            void set(const char *s, sindex *namei = nullptr);
            String toString(stringarray *names = nullptr) const;
            void init();
            
            bool operator < (const sbpos &p) const;
            bool operator == (const sbpos &p) const;
            bool operator != (const sbpos &p) const;
        };
        /*
        struct sbexpos : public sbpos {
            String posid;
            float score;
            sattribute attribute;
        };
        
        class SBPosList : public Array<sbexpos> {
        private:
            
        public:
            
        };
        */
    }
}

#endif
