#ifndef __SBIO_SVARIANT_H__
#define __SBIO_SVARIANT_H__

#include "sutil/sthread.h"
#include "sbioinfo/sbannot.h"
#include "sbioinfo/sbioseq.h"
#include "sbioinfo/svarutil.h"

namespace slib {
    namespace sbio {
        
        #define SMALL_VARIANT 0x01
        #define CN_VARIANT 0x02
        #define SR_VARIANT 0x04
        
#define DEFAULT_BG_DEPTH 3.0
#define DEFAULT_DELCP 0.75
#define DEFAULT_DUPCP 1.5
#define DEFAULT_MULCP 2.5
#define DEFAULT_HOMO_DEL 0.35
#define DEFAULT_HOMO_FREQ 0.65
#define DEFAULT_CP_DIFF 0.25
#define DEFAULT_MIN_CVLEN 20
#define DEFAULT_CP_PROB 0.01
        
#define DEFAULT_SR_READ 3
#define DEFAULT_MIN_QUAL 1.0
#define DEFAULT_MIN_FREQ 0.01
#define DEFAULT_FR_BIAS 1.0
#define DEFAULT_COMP_BIAS 1.0
        
#define DEFAULT_VC_DP 5
#define DEFAULT_VC_VDP 3
#define DEFAULT_VC_FREQ 0.1
#define DEFAULT_VC_HOMO_FREQ 0.75
#define DEFAULT_VC_QUAL 10.0
        
#define DEFAULT_MIN_VAR_SIZE 1
#define DEFAULT_MAX_VAR_SIZE -1
#define DEFAULT_DIST 50
#define DEFAULT_DIFF 100
#define DEFAULT_FREQ_BIN 20
        
        struct svc_param {
            //0:SNV, 1:MNV, 2:DEL, 3:INS
            sint min_depth[4], min_vdepth[4];
            double min_freq[4], homo_freq[4], min_qual[4];
            bool repeat_base, homo_select,
            cds_var, utr_var, exon_var, intron_var, spsite_var;
            
            svc_param();
            svc_param(const svc_param &p);
            ~svc_param();
            
            svc_param &operator=(const svc_param &par);
            void set(const sobj &obj);
            sobj toObj();
        };
        
        struct scnv_param {
            float bg_depth, del_cp, dup_cp, mul_cp, homo_del_cp, homo_freq, cp_diff, min_cv_len;
            
            scnv_param();
            scnv_param(const scnv_param &p);
            ~scnv_param();
            
            scnv_param &operator=(const scnv_param &par);
            void set(const sobj &obj);
            sobj toObj();
        };
        
        struct ssrv_param {
            bool detect_var[5], detect_comp_var[4];
            sint min_sr[5], min_comp_sr[4];
            double min_freq, max_fr_bias, max_comp_bias;
            
            ssrv_param();
            ssrv_param(const ssrv_param &p);
            ~ssrv_param();
            
            ssrv_param &operator=(const ssrv_param &par);
            void set(const sobj &obj);
            sobj toObj();
        };
        
        struct svariant_param {
            sint min_length[3], max_dist, max_diff, min_qual;
            
            scnv_param cnv_par;
            ssrv_param srv_par;
            svc_param smv_par;
            
            svariant_param();
            svariant_param(const svariant_param &p);
            ~svariant_param();
            
            void set(const sobj &obj);
            sobj toObj();
        };
        
        
        extern double read_bias(sint r1, sint r2);
        extern double phred_val(double q);
        
        class SVariant;
        class SVarList;
        class svcf;
        
        struct svar_data {
            sushort type;
            sbpos pos1, pos2;
            String alt;
            sushort pread, nread;
            double qual;
            
            svar_data();
            svar_data(const sbpos &pos);
            svar_data(const svar_data &v);
            ~svar_data();
            
            svar_data &operator = (const sbpos &pos);
            svar_data &operator = (const svar_data &v);
            svar_data &operator += (const svar_data &v);
            
            svar_data &dat();
            const svar_data &dat() const;
            
            void comp();
            sint total() const;
            double bias() const;
            double phred() const;
            bool comparable(const svar_data &v) const;
            bool lt(const svar_data &var, size_t dist) const;
            bool equal(const svar_data &var, size_t dist) const;
            
            bool operator < (const svar_data &v) const;
            bool operator ==(const svar_data &v) const;
        };
        
        typedef Array<svar_data> vararray;
        
        /*
         small var ... consequence
         large var ... null, dup/mul, fusion
         
         conseq... amino acid change, frame shift/in frame per transcript
         null/dup
         fusion... hetero(ipsi/contra tandem/internal), homo(tandem, inv, del)
         
         *gene_id, gene_name, gene_dir, mutated site, mut_type
         aachange {ori alt pos}, fusion {target type}
         
         type:
         substitute
            Silent  0
            Mis     1
            Non     2
         sizevar
            FrameShift   4
            InFrame      8
         
            Null    10
            Dup/Mul 20
         
         Rearrange  40
         Fusion     80
         
         tandem      1
         ins         2
         wrap        4
         ipsi        10
         contra      20
         
         ->...> tandem ipsi
         -> <... tandem contra
         -> ...> -> wrap/ins ipsi
         -> <... -> wrap/ins contra
         
         */
        
        class SVariant : public svar_data {
			friend SVarIO;
            friend SVarList;
			friend SVarAnnotator;
			friend SVarFilter;
            
        public:
            struct transcript_site {
                sushort type, site;
                String trs_name, ori, alt;
                sint pos;
                
                transcript_site();
                transcript_site(sbio::SBAnnotDB::transcript_info *ti);
                transcript_site(const SVariant::transcript_site &trs);
                ~transcript_site();
                
                transcript_site &operator=(const SVariant::transcript_site &trs);
            };
            
            struct gene_site {
                sushort type;
                String gene_name;
                bool gene_dir;
                slib::Array<transcript_site> transcripts;
                //intarray fusion;
                
                gene_site();
                gene_site(const sbio::SBAnnotDB::gene_info *gi);
                gene_site(const SVariant::gene_site &g);
                ~gene_site();
                
                gene_site &operator=(const SVariant::gene_site &g);
            };
            
        public:
            static String vmethod(subyte i);
            static subyte vmethodIdx(const char *s);
            
            static String vsite(sushort i);
            static subyte vsiteIdx(const char *s);
            static String mtype(sushort i);
            static subyte mtypeIdx(const char *s);
            static String vannot(Array<SVariant::gene_site> &genes);
            
            
        public:
			sushort flag;
            subyte method;
            String name, ref_id1, ref_id2;
            double copy[6], bg_copy[4], frequency;
            bool homo;
            slib::Array<SVariant::gene_site> genes;
            stringarray mutants;
            SDictionary attribute;
            
        public:
            SVariant();
            SVariant(subyte m, const svar_data &v);
            SVariant(const SVariant &var);
            ~SVariant();
            
            SVariant &operator=(const SVariant &var);
            
            bool comparable(const SVariant *var) const;
            bool lt(const SVariant *var, size_t dist) const;
            bool equal(const SVariant *var, size_t dist) const;
            
            bool operator <(const SVariant &var) const;
            bool operator ==(const SVariant &var) const;
        };
        
#define svar sptr<SVariant>
        
        class SVarList : public Array<svar> {
			friend SVarIO;
			friend SVariant;
			friend SVarAnnotator;
			friend SVarFilter;

            typedef std::pair<sbyte, const SVariant *> ivpair;
            typedef RArray<ivpair> ivarray;
            
        public:
            static const sattribute ION_VAR_CALL;
            
        protected:
            svariant_param *par;
            sindex ref_name_idx;
            size_t map_error_dist;
            
        public:
            String vl_name;
            SDictionary vl_attribute;
            
        private:
            void _annotate(SVariant *var, SBSeqList *ref, SBAnnotDB *annot);
            void _shrink(size_t s);
            void _shrink(ivarray &array, size_t s);
            
        public:
            SVarList(const char *s = nullptr);
            SVarList(const SVarList &vl);
            virtual ~SVarList();
            
            void loadVCF(const char *path, SBSeqList *ref, const sattribute &converter = {});
            
            void load(const char *path);
            
            void saveVCF(const char *path, SBSeqList *ref);
            void saveTable(const char *path);
            void save(const char *path);
            
            void setParam(svariant_param *p);
            void setReference(sbio::SBSeqList *r);
            
            void removeDup();
            void merge(const SVarList &vl);
            void subtract(const SVarList &vl);
            void unique(SVarList &uni, const SVarList &vl);
            void common(SVarList &com, const SVarList &vl);
            void annotate(SBAnnotDB *annot, sbio::SBSeqList *ref, stpool *threads = nullptr);
            
			void tidy(size_t &s);

            sindex &nameIndex();
            
            void clearAll();
        };

		
        
        class svcf : public SVarList {
        private:
            bool format;
            sattribute *converter;
            
        public:
            svcf();
            ~svcf();
            
            void importVCF(const char *path, sattribute *converter = nullptr);
            void exportVCF(const char *path);
            void setConverter(sattribute *c);
            
            void checkRepeat();
        };
    }
}
#endif
