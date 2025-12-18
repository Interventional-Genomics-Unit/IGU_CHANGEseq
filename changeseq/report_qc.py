import os
import logging


logger = logging.getLogger('root')
logger.propagate = False


def write_qc(qc_file,preprocessed_logfile,coverage_stat_file):
    #for sample in c.parameters['samples']:
    # preprocessed_logfile = os.path.join(c.parameters["analysis_folder"], 'preprocessed', sample + '.txt')
    # coverage_stat_file  = os.path.join(c.parameters["analysis_folder"], 'qc', sample + '_aligned_stats.txt')
    #
    # qc_file = os.path.join(c.parameters["analysis_folder"], 'qc', sample + '_qc_report.txt')
    metrics = ['total_reads',
               'reads_tn5_filtered',
               'preprocessed_read_remaining',
               'reads_aligned',
               'reads_mapped_outie - derived from circles',
               'total_reads_mapped_inner - derived from linear',
               'mate mapped to different chromosome',
               'average read length',
               'average mapping quality']
    metrics_dict =dict()
    metrics_dict.fromkeys(metrics)
    if os.path.isfile(preprocessed_logfile):
        lines = open(preprocessed_logfile,'r').readlines()


        for i,line in enumerate(lines[0:100]):
            if 'Total read pairs processed:' in line:
                total = int(line.split(':')[-1].strip().replace(",", ""))*2
                metrics_dict['total_reads'] = total

            if 'Read 1 with adapter:' in line:
                r1_n,r1_p = line.split(':')[-1].strip().split(" ")
                r2_n,r2_p = lines[i+1].split(':')[-1].strip().split(" ")
                n = int(r1_n.replace(",", "")) +  int(r2_n.replace(",", ""))
                metrics_dict['reads_tn5_filtered']= f"{n} {round((n/total)*100,1)}%"


                n,percent = lines[i+3].split(':')[-1].strip().replace(",", "").split(" ")
                percent =percent.replace(")", "").replace("(", "")
                n =int(n)*2
                metrics_dict['preprocessed_read_remaining']= f"{n} {percent}"

        if os.path.isfile(coverage_stat_file):
            lines = open(coverage_stat_file,'r').readlines()

            for i,line in enumerate(lines[0:100]):
                if 'reads mapped:' in line:

                    n_aligned = int(line.split(':')[-1].strip())
                    metrics_dict[ 'reads_aligned'] = f"{n_aligned} {round((n_aligned/total)*100,1)}%"
                if 'outward oriented pairs' in line:
                    out_aligned = int(line.split(':')[-1].strip())
                    metrics_dict['reads_mapped_outie - derived from circles'] =  f"{out_aligned} {round((out_aligned/total)*100,1)}%"
                if 'inward oriented pairs' in line:
                    in_aligned = int(line.split(':')[-1].strip())
                    metrics_dict['total_reads_mapped_inner - derived from linear'] = f"{in_aligned} {round((in_aligned/total)*100,1)}%"
                if 'pairs on different chromosomes:' in line:
                    n = int(line.split(':')[-1].strip())
                    metrics_dict['mate mapped to different chromosome'] = f"{n} {round((n/total)*100,1)}%"
                if 'average length' in line:
                    n = line.split(':')[-1].strip()
                    metrics_dict['average read length'] = n
                if 'average quality' in line:
                    n = line.split(':')[-1].strip()
                    metrics_dict['average mapping quality'] = n


        with open(qc_file, "w") as out:
            for k,v in metrics_dict.items():
                out.write(f"{k}: {v}\n")

