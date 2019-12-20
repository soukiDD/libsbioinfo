#ifndef SBIO_SBIOANNOT_H
#define SBIO_SBIOANNOT_H

#include "sdb/sdatabase.h"
#include "sbioinfo/sbioutil.h"

namespace slib {
    using namespace sio;
    
    namespace sbio {
        class SBAnnotDB : public sdb::SDataBase {
            /*
             SQLite Table:
                INFO:
                    ID(KEY)
                    NAME(STRING)
                        e.g. creator
                             date
                             version
                             species
                             etc...
             
                    VALUE(STRING)
             
                CHROMOSOME:
                    ID(KEY)
                    NAME(STRING)
                    LENGTH(INTEGER)
                
                CONTIG:
                    ID(KEY)
                    NAME(STRING)
                    CHROMOSOME(INTEGER)
                    START(INTEGER)
                    END(INTEGER)
             
                GENE:
                    ID(KEY)
                    GENE_ID(STRING)
                    TYPE(INTEGER)
                    NAME(STRING)
                    CHROMOSOME(INTEGER)
                    START(INTEGER)
                    END(INTEGER)
                    STRAND(STRING)
                    OTHER_NAME(STRING)
                    TRANSCRIPT(ARRAY)
                    DESCRIPTION(STRING)
                    ATTRIBUTE(STRING)
                
                TRANSCRIPT:
                    ID(KEY)
                    GENE_ID(INTEGER)
                    TYPE(INTEGER)
                    CHROMOSOME(INTEGER)
                    NAME(STRING)
                    START(INTEGER)
                    END(INTEGER)
                    STRUCTURE(ARRAY)
             
                STRUCTURE:
                    ID(KEY)
                    TRANSCRIPT_ID(INTEGER)
                    TYPE(INTEGER)
                    START(INTEGER)
                    END(INTEGER)
                    
                MUTATION:
                    ID(KEY)
                    MUT_ID(STRING)
                    TYPE(INTEGER)
                    NAME(STRING)
                    CHROMOSOME(INTEGER)
                    START(INTEGER)
                    END(INTEGER)
                    STRAIN(STRING)
                    ATTRIBUTE(STRING)
             
                VARIATION:
                    ID(KEY)
                    VAR_ID(STRING)
                    TYPE(INTEGER)
                    NAME(STRING)
                    CHROMOSOME(INTEGER)
                    START(INTEGER)
                    END(INTEGER)
                    STRAIN(STRING)
                    ATTRIBUTE(STRING)

                FEATURE:
                    ID(KEY)
                    TYPE(INTEGER)
                    NAME(STRING)
                    START(INTEGER)
                    END(INTEGER)
                    ATTRIBUTE(ARRAY)
             */
            
        public:
            struct annot_info : public sbpos {
                sint _id, type;
                String name;
                
                annot_info();
                annot_info(const sint &t, const String &n, const sbpos &p);
                annot_info(const annot_info &info);
                ~annot_info();
            };
            
            struct chr_info : public annot_info {
                chr_info();
                chr_info(const char *s, const sbpos &pos);
                chr_info(const chr_info &c);
                ~chr_info();
                chr_info &operator=(const chr_info &c);
            };
            
            struct contig_info : public annot_info {
                contig_info();
                contig_info(const contig_info &c);
                ~contig_info();
                contig_info &operator=(const contig_info &c);
            };
            
            struct struct_info : public annot_info {
                struct_info();
                struct_info(const struct_info &s);
                ~struct_info();
                
                struct_info &operator=(const struct_info &s);
            };
            struct gene_info;
            struct transcript_info : public annot_info {
                gene_info *gene;
                Array<struct_info> structures;
                
                transcript_info();
                transcript_info(const transcript_info &t);
                ~transcript_info();
                
                transcript_info &operator=(const transcript_info &t);
                
                void addStructure(struct_info &&s);
                void setGene(gene_info *g);
                
                sregion exonRegion();
                sregion codingRegion();
                Array<struct_info> messenger();
            };
            
            struct gene_info : public annot_info {
                String gene_id, description;
                stringarray other_names;
                CArray<transcript_info *> transcripts;
                sattribute attribute;
                
                gene_info();
                gene_info(const gene_info &g);
                ~gene_info();
                
                gene_info &operator=(const gene_info &g);
                
                void setDescription(String *str);
                void addTranscript(transcript_info *t);
            };
            
            struct variant_info : public annot_info {
                String var_id, strain;
                sattribute attribute;
                
                variant_info();
                variant_info(const variant_info &v);
                ~variant_info();
                
                variant_info &operator=(const variant_info &v);
            };
            
            struct feature_info : public annot_info {
                feature_info();
                feature_info(const feature_info &f);
                ~feature_info();
                
                feature_info &operator=(const feature_info &f);
            };
            
            static const int16_t LOAD_CHR = 0x01;
            static const int16_t LOAD_CTG = 0x02;
            static const int16_t LOAD_GENE = 0x04;
            static const int16_t LOAD_TRANS = 0x08;
            static const int16_t LOAD_MUT = 0x10;
            static const int16_t LOAD_VAR = 0x20;
            static const int16_t LOAD_FTR = 0x40;
            
            typedef Array<SBAnnotDB ::chr_info> chr_vec;
            typedef Array<SBAnnotDB ::contig_info> ctg_vec;
            typedef CArray<SBAnnotDB ::contig_info *> ctgp_vec;
            typedef Array<SBAnnotDB ::gene_info> gene_vec;
            typedef CArray<SBAnnotDB ::gene_info *> genep_vec;
            typedef Array<SBAnnotDB ::transcript_info> trs_vec;
            typedef CArray<SBAnnotDB ::transcript_info *> trsp_vec;
            typedef Array<SBAnnotDB ::variant_info> var_vec;
            typedef CArray<SBAnnotDB ::variant_info *> varp_vec;
            typedef Array<SBAnnotDB ::feature_info> ftr_vec;
            typedef CArray<SBAnnotDB ::feature_info *> ftrp_vec;
            
            typedef std::pair<String *, sint> name_pair;
            
        private:
            int16_t _mode;
            Array<sorder> bin_order;
            sindex _chr_index;
            Array<name_pair> _gene_name_index;
            Array<Array<ctgp_vec>> _ctg_index;
            Array<Array<genep_vec>> _gene_index;
            Array<Array<trsp_vec>> _trs_index;
            Array<Array<varp_vec>> _mut_index, _var_index;
            Array<Array<ftrp_vec>> _ftr_index;
            
            intarray _ctg_name_index, _trs_name_index, _mut_name_index, _var_name_index, _ftr_name_index;
            
        public:
            chr_vec chromosomes;
            ctg_vec contigs;
            gene_vec genes;
            trs_vec transcripts;
            var_vec mutants, variations;
            ftr_vec features;
            
        private:
            void initIdx();
            void loadChrInfo();
            void loadContigInfo();
            void loadGeneInfo();
            void loadTranscriptInfo();
            void loadMutantInfo();
            void loadVariationInfo();
            void loadFeatureInfo();
            
        public:
            SBAnnotDB ();
            SBAnnotDB (const char *path, int m = LOAD_CHR);
            ~SBAnnotDB ();
            
            void open(const char *path);
            void load(int m = LOAD_CHR);
            int mode() const;
            void setMode(int m);
            
            String species();
            String version();
            
            size_t chrNum() const;
            SBAnnotDB ::chr_info chrInfo(int idx) const;
            size_t chrIndex(const char *name) const;
            const sindex &nameIdx() const;
            
            void ctgInfo(ctgp_vec &vec, const sbpos &pos);
            void ctgInfo(ctgp_vec &vec, const char *name, sdb::SQL_MATCH match = sdb::EXACT);
            
            void geneInfo(genep_vec &vec, const sbpos &pos, bool trans = true);
            void geneInfo(genep_vec &vec, const char *name, bool trans = true, sdb::SQL_MATCH match = sdb::EXACT);
            void nearestGeneInfo(genep_vec &vec, const sbpos &pos, bool trans = true);
            
            void transcriptInfo(trsp_vec &vec, const sbpos &pos, bool gene = false);
            void transcriptInfo(trsp_vec &vec, const char *name, bool gene = false, sdb::SQL_MATCH match = sdb::EXACT);
            
            void mutantInfo(varp_vec &vec, const sbpos &pos);
            void mutantInfo(varp_vec &vec, const char *name, sdb::SQL_MATCH match = sdb::EXACT);
            
            void variationInfo(varp_vec &vec, const sbpos &pos);
            void variationInfo(varp_vec &vec, const char *name, sdb::SQL_MATCH match = sdb::EXACT);
            
            void featureInfo(ftrp_vec &vec, const sbpos &pos);
            void featureInfo(ftrp_vec &vec, const char *name, sdb::SQL_MATCH match = sdb::EXACT);
        };
    }
}
#endif
