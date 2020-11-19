# python version of control pipelines script 
# Few comments: 
# - script can submit many many jobs. You need to monitor it. 
# - Currently doesn't take into account the mysql servers load - but the analyses that load the servers. 
# - script consider number of connections to mysql servers. 

#!/usr/bin/python




# pip install colorama
import sys
import argparse
import logging
import subprocess
import sqlalchemy as db
from subprocess import Popen, PIPE
from sqlalchemy.sql import select
from sqlalchemy import func, and_
from urllib.parse import urlparse
from sqlalchemy.orm import sessionmaker
from colorama import init
from colorama import Fore, Back, Style

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql
from numpy.f2py.crackfortran import verbose
from test.test_pprint import set2
pymysql.install_as_MySQLdb()




class Pipeliner:

    # Constructor
    def __init__(self, name, pipeline_status, mysql_status=0, farm_status=0, pipeline_input_ids_status=0):
        self.name            = name  
        self.pipeline_status = pipeline_status
        colour_pipeline(self.name, self.pipeline_status)
        self.mysql_status    = mysql_servers_status(self.name, 2, True)
        self.farm_status     = check_farm_jobs() 
        self.pipeline_input_ids_status = check_pipeline_inputIDs_status(self.name, verbose)

    def print_data(self, verbose=False):
        if verbose:
            print('# name, {}'.format(self.name))
            print('# pipeline status, {}'.format(self.pipeline_status))
            print('# mysql_status, {}'.format(self.mysql_status))
            print('# Can I submit more?, {}'.format(self.farm_status))
            print('# pipeline_input_ids_status, {}'.format(self.pipeline_input_ids_status))
    
    def calculate_bkeeper_params(self, number_of_loops, sleep_div, number_of_workers, verbose = True): 
        max_of_mysql_connections = sum(self.mysql_status.values())
        
        verbose=True
        if verbose: 
            print('# script will do the calculations')
            print('max_of_mysql_connections:', max_of_mysql_connections)

        ### MAJOR CHECKS ### 
        # STOP/GO checks
        while max_of_mysql_connections  > 2000:
            print('max_connections are: very high (', max_of_mysql_connections , ')', '. You need to wait.' )
            sleep(10)
            max_of_mysql_connections  = sum(mysql_servers_status(self.name, check=2))
        
        while not self.farm_status:
            print('total_number_of running jobs (in farm) are: Very high. You need to wait.' )
            sleep(10)
            self.farm_status = check_farm_jobs() 
            # keep sending jobs else sleep.     

        if  ( ( self.pipeline_status == 'boost_run') and ( self.pipeline_input_ids_status['READY'] > 1000 ) ): 
            # this is going to be a real boost for parts of high priority pipelines that makes connections to mysql servers.  
            print('I am in a boost mode here. I will take everything available dude')
            max_loops = 15
            sleep_div = 0.2            
            how_much_space = 2500 - max_of_mysql_connections
            number_of_workers = number_of_workers + (how_much_space/2) # how_much_space/2 is the high-priority bonus workers
            print('I have: ', how_much_space, ' to use and I will use: ',  number_of_workers)
        elif ( (max_of_mysql_connections < 500 ) and (self.pipeline_input_ids_status['READY'] > 10000 ) ):
            print('+10000_ready_jobs_are_waiting')
            number_of_loops = 15
            sleep_div = 0.2
            number_of_workers = number_of_workers*1.8
        elif ( (max_of_mysql_connections < 1000 ) and ( self.pipeline_input_ids_status['READY'] > 10000 ) ):
            number_of_loops = 8
            sleep_div = 0.3
            number_of_workers = number_of_workers*1.5
        elif ( (max_of_mysql_connections < 1000 ) and ( self.pipeline_input_ids_status['READY'] > 1000 ) ): 
            print('normal')
            number_of_loops = 5
            sleep_div = 0.1 
            number_of_workers = number_of_workers*1.2            
        elif ( (max_of_mysql_connections < 1000 ) and ( self.pipeline_input_ids_status['READY'] < 5 ) 
               and ( self.pipeline_input_ids_status['READY'] > 0 ) ):
            print('not much to run, I will loop twice')
            number_of_loops = 2
        elif ( (max_of_mysql_connections < 1000 ) and ( self.pipeline_input_ids_status['READY']) ):
            number_of_loops = 10
        elif ( self.pipeline_input_ids_status['READY'] == 0 ) :
            if ( self.pipeline_input_ids_status['RUN'] == 0 ) : 
                print('Nothing_to_run. Looks like we might have finish! ')
                number_of_loops = 1                
        else: 
            print('Problem: Something else is going on. ')
        # number_of_loops = 1
        
        
        return number_of_loops, sleep_div, number_of_workers

    
def ready_analyses(db_url, verbose=False):
    engine = db.create_engine(db_url)
    connection = engine.connect()
    connections_info = {}
    
    if connection:
        metadata = db.MetaData()
        j = db.Table('job',metadata, autoload=True, autoload_with=engine)
        ab = db.Table('analysis_base',metadata,autoload=True,autoload_with=engine)
        s = select([j.columns.analysis_id, ab.columns.logic_name]).\
        where(and_((j.c.analysis_id == ab.c.analysis_id),(j.c.status == 'READY')))
        results = connection.execute(s)
        
        for analysis_id, logic_name in results: 
            connections_info[logic_name] = analysis_id
            print('ready analysis:',analysis_id, 'logic_name', logic_name)
    else:
        raise ValueError(f'Could not connect to the target database {db_url}.')
    return connections_info
    

# run commands 
def block_checkcall(cmd, verbose) :
    if verbose:
        print('running cmd: ', cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    if verbose == 3: 
        print("err:", stderr)
        print("out:", stdout)
    return stdout, stderr 


# check run/pending/fail jobs: 
def check_farm_jobs(task_max=2000, verbose=True):
    stdout_run, stderr_run = block_checkcall('bjobs | grep RUN | wc -l', verbose)
    njobs_run = int(stdout_run.rstrip())
    stdout_pend, stderr_pend = block_checkcall('bjobs | grep PEND | wc -l', verbose)
    njobs_pend = int(stdout_pend.rstrip())
    if verbose: 
        print("pending_jobs:",njobs_pend)
        print("running_jobs:",njobs_run)
        
    if njobs_run > task_max:
        return False
    else: 
        print('You are running less than ', task_max , ' jobs')
        return True

def mysql_servers_status(db_url, check=1, verbose=False):
    url_attributes = urlparse(db_url)
    connections_info = {}

    if check == 1:
        list_cmd = [url_attributes.hostname , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
        cmd_processlist = ''.join(list_cmd)
        stdout_pend, stderr_pend = block_checkcall( cmd_processlist , verbose)
        connections_info['total'] = stdout_pend
    elif check == 2:
        # This will calculate average
        set1 = ('mysql-ens-genebuild-prod-2','mysql-ens-genebuild-prod-3', 'mysql-ens-genebuild-prod-4')
        set2 = ('mysql-ens-genebuild-prod-5','mysql-ens-genebuild-prod-6', 'mysql-ens-genebuild-prod-7')
        max_connections = int(1)
        if url_attributes.hostname == 'mysql-ens-genebuild-prod-4' : 
            for test_host in set1:
                list_cmd = [test_host, "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
                cmd_processlist = ''.join(list_cmd)
                stdout_pend, stderr_pend = block_checkcall( cmd_processlist , verbose)
                stdout_pend = int(stdout_pend)
                if max_connections < stdout_pend: 
                    max_connections = stdout_pend
        elif url_attributes.hostname == 'mysql-ens-genebuild-prod-7' :
            for test_host in set2:
                list_cmd = [test_host , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
                cmd_processlist = ''.join(list_cmd)
                stdout_pend, stderr_pend = block_checkcall( cmd_processlist , verbose)
                stdout_pend = int(stdout_pend)
                if max_connections < int(stdout_pend): 
                    max_connections = stdout_pend                
        else :
            print('problem')
        
        connections_info['max'] = max_connections

    else: 
        list_cmd = [ url_attributes.hostname , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | ", 
                    " awk ", " '{print $5}' "  ]
        cmd_processlist = ''.join(list_cmd)
        stdout_pend, stderr_pend = block_checkcall( cmd_processlist , verbose)

        for ser_status in stdout_pend.split(): 
            if (ser_status in connections_info): 
                connections_info[ser_status] += 1 
            else: 
                connections_info[ser_status] = 1 

    if verbose: 
        for xx in connections_info:
            print (xx,':',connections_info[xx])
       
        
    return connections_info 
    

def check_pipeline_inputIDs_status(db_url, verbose = False):
    freq = {}
    engine = db.create_engine(db_url)
    connection = engine.connect()
    if connection:
        metadata = db.MetaData()
        job = db.Table('job',metadata,autoload=True,autoload_with=engine)
        query = db.select([job.columns.status,job.columns.job_id])
        query = query.select_from(job)
        results = connection.execute(query).fetchall()
        for status,job_id in results: 
            if (status in freq): 
                freq[status] += 1 
            else: 
                freq[status] = 1 
    
        if verbose: 
            for xx in freq:
                print (xx,':',freq[xx])
    else:
        raise ValueError(f'Could not connect to the target database {db_url}.')
    return freq
  
  
def colour_pipeline(pipeline_name, pipeline_stage): 
    init(autoreset=True)
    if (pipeline_stage in ['run', 'boost_run'] ): 
        print(Fore.GREEN + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'stop'):
        print(Fore.MAGENTA + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'failed'):    
        print(Fore.RED + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'ready'):
        print(Fore.YELLOW + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'done'):
        print(Fore.WHITE + pipeline_name , pipeline_stage )
    else: 
        print('are you sure about the status?') 

        
# check boost priority
def use_priority(pipelines):
    boost_option = input('Do you want to boost pipelines(y/n)?')
    if boost_option.startswith('y'): 
        for pipeline_priority in pipelines:
            boost_pipeline = input('boost pipeline:' + pipeline_priority + 'y/n?') 
            if boost_pipeline.startswith('y'):
                pipelines[pipeline_priority] = 'boost_run'
     
    return pipelines            

def reset_priorities(pipelines): 
    print('Enough priorities')
    for pipeline_priority in pipelines:
        boost_pipeline = pipelines[pipeline_priority]
        if boost_pipeline.startswith('boost_run'):
            pipelines[pipeline_priority] = 'run'
    
    return pipelines

# Analyses that don't push mysql servers a lot and need to submit many jobs. 
# It is not worth to put here an analysis that submits only 1-2 jobs      
def mysql_load_analyses(list_type):
    if list_type == 'mysql_Free': 
        analyses_mysql_free = ('download_long_read_fastq', 'minimap2', 'minimap2_himem', 'run_repeatmasker',\
                     'rebatch_repeatmasker','run_repeatmasker_small_batch','download_RNASeq_fastqs','bwa',\
                     'bwa2bam','merged_tissue_file','bam2bigwig')
        return analyses_mysql_free
    elif list_type == 'mysql_Load': 
        # analyses that load mysql servers (experience speaking)
        analyses_mysql_load = ('layer_annotation','run_utr_addition','genebuilder')
        return analyses_mysql_load
    elif list_type == 'quickly_done':
        # analyses that are quickly done and can be loop more times when in high priority
        # https://docs.google.com/spreadsheets/d/1RgPW9q05YGENajVs9XsEr_XN7qZCM7rGt-_d0z3iZ1Y/edit?usp=sharing
        analyses_that_run_quickly = ('core_assembly_name_update', 'otherfeatures_assembly_name_update', 'pairaligner_stats',\
                                      'get_species_list', 'add_method_link_species_link_tag', 'rebatch_repeatmasker',\
                                      'update_rnaseq_ise_logic_names', 'load_selenocysteine', 'create_merge_analyses_type_job',\
                                      'fan_rnaseq_for_layer_db', 'fan_projection', 'null_otherfeatures_columns',\
                                      'create_tissue_jobs', 'create_overlapping_slices', 'email_loading_report',\
                                      'final_meta_updates', 'fan_genscan_blasts', 'fan_genblast_rnaseq_support',\
                                      'semaphore_10mb_slices', 'create_toplevel_input_ids', 'generate_besttargetted_jobs',\
                                      'update_cdna_analyses', 'create_cdna2genome_slices', 'cesar_failed_projection',\
                                      'failed_ig_tr_genblast_proteins', 'create_repeatmasker_slices', 'drop_backup_tables_job',\
                                      'create_analyses_type_job', 'clean_unused_analyses', 'update_ig_tr_biotypes',\
                                      'split_ig_tr_genblast_jobs', 'update_ig_tr_hitnames', 'fan_refseq_import', 'fan_lastz',\
                                      'indicate_proteome', 'indicate_BT', 'create_genscan_slices', 'load_meta_info',\
                                      'create_cdna_toplevel_slices', 'backup_original_csv', 'fan_long_read',\
                                      'create_genblast_rnaseq_slice_ids', 'create_seleno_homology_jobs', 'checking_file_path',\
                                      'fan_ncrna', 'create_fastq_download_jobs', 'remove_rnaseq_for_layer_daf_features',\
                                      'set_otherfeatures_meta_levels', 'fan_merge_analyses', 'create_header_intron',\
                                      'dump_features', 'generate_ig_tr_jobs', 'restore_ig_tr_biotypes', 'generate_pmatch_jobs',\
                                      'generate_targetted_jobs', 'update_otherfeatures_db', 'create_10mb_slice_ids',\
                                      'change_biotype_for_weak_cds')
        return analyses_that_run_quickly
    else:
        print('problem')    


# connect to registry db and update status of the pipeline





def main():
    """
    main function
    """
    # get arguments from command line: 
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--input", help="Your input file.")
    parser.add_argument("-a", "--act", help="What I have to do?")
    parser.add_argument("-l", "--number_of_loops", type=int, help="How many loops?", default=10 )
    parser.add_argument("-w", "--number_of_total_workers", type=int, help="How many total_workers?", default=250 )
    parser.add_argument("-s", "--skip_checks", help="Skip tests?", default=False, action='store_true')
    parser.add_argument("-v", "--verbose",dest='verbose',action='store_true', help="Verbose mode.")
    options = parser.parse_args()
    verbose = options.verbose
    if verbose:
        print("Verbose mode on")
    else:
        print("Verbose mode off")

    # standard parameters: 
    beekeeper_locations = '/nfs/production/panda/ensembl/kbillis/enscode_2020_08/enscode/ensembl-hive/scripts/beekeeper.pl '
    debug_parameter     = ' -debug 1 '
    how_many_loops      = options.number_of_loops
    priority_limit      = int(how_many_loops*0.1)
    skip_checks         = options.skip_checks


    ### Start running the script: ### 
    # read a dictionary with pipelines to run: 
    d = {}
    with open(options.input) as f:
        data = f.readlines()
        for line in data:
            if line.startswith('#'):
                print('WARN: You comment out this pipeline: ',line)
            else:
            # key = line.split()
                print(line)
                d[line.rstrip()] = 'run'
    
    
    if options.act == 'loop':
        d = use_priority(d)
    else: 
        how_many_loops = 1
        
    # how many total loops 
    # this is a while: 
    for ii in range(how_many_loops): 
        # if ii == priority_limit:
        #    d = reset_priorities(d)
        
        # check how many workers each pipeline should run. 
        bkeeper_total_workers = options.number_of_total_workers
        worker_equal_dis =  bkeeper_total_workers/len(d)
        
        
        ## loop those pipelines 
        for pipeline_name in d:
            bkeeper_max_loops = 5
            bkeeper_sleep = 0.5

            print('\n#################################################################')
            print(pipeline_name,'status',':',d[pipeline_name])
            
            if options.act != 'loop':
                mylist = ['perl ', beekeeper_locations, ' ', debug_parameter , ' -' , options.act , ' -url ' , pipeline_name ]
                cmd_tmp = ''.join(mylist)
                stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose)
                # Construct an instance of the Pipeline class 
                g = Pipeliner(pipeline_name, d[pipeline_name])  
                g.print_data(verbose=True)   # Call an instance method; prints 
            else: 
                if ( d[pipeline_name] == 'failed') :
                    print('This pipeline has FAILED issues:', d[pipeline_name])
                    colour_pipeline(pipeline_name, d[pipeline_name]) 
                    continue
                elif ( d[pipeline_name] == 'failed' ) :
                    print('I will create config')
                    # method to create config
                    d[pipeline_name] = 'run'
                    continue
                
                analyses_ready = ready_analyses(pipeline_name)
                analysis_list = mysql_load_analyses(list_type='mysql_Free')
                is_mysql_Free = 'no'
                for logic_name in analyses_ready: 
                    if logic_name in analysis_list: 
                        print('yes!', logic_name)
                        is_mysql_Free = 'yes' 
                if skip_checks: 
                     bkeeper_total_workers = worker_equal_dis
                elif is_mysql_Free == 'yes': 
                     bkeeper_total_workers = 500
                     bkeeper_max_loops = 10
                else:
                    # Here is the semi-dynamic part
                    # Construct an instance of the Pipeline class
                    print('I will set beekeeper: sleep, max_loops and total_running_workers_max')
                    g = Pipeliner(pipeline_name, d[pipeline_name])
                    # g.print_data(verbose=True)  
                    # will return number of workers, amount of loops and sleep value.
                    bkeeper_max_loops, bkeeper_sleep, bkeeper_total_workers = g.calculate_bkeeper_params(bkeeper_max_loops, bkeeper_sleep, worker_equal_dis)
                    
                mylist = ['perl ' , beekeeper_locations , ' ' , debug_parameter , ' -' , options.act , ' -url ' , pipeline_name ,\
                          ' -max_loops ' , str(bkeeper_max_loops) , ' -total_running_workers_max ' , str(bkeeper_total_workers) ,\
                          ' -sleep ' , str(bkeeper_sleep) ]
                cmd_tmp = ''.join(mylist)
                print('command: ',cmd_tmp)
                # STOP/GO checks should be here too.  
                # stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose)
                
                if (d[pipeline_name] == 'boost_run'):
                    while check_pipeline_inputIDs_status(pipeline_name)['READY'] > 5000: 
                        # give it a bit more of life until more input_id are done.
                        # stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose) 
                        print('I like this pipeline, ', pipeline_name, 'and I will run again. ')
                        sleep(5)
            
            # TODO: major check if there is an analysis that takes too much time (for example, if one analysis runs for more than 10 hours,
            # mark db as 'problematic/failed'): 
            

             
        
        # call the class 
        # pipeline, mysql_check, check_jobs 
    print('Game over!')


if __name__ == "__main__":
    main()

