# Python version of control pipelines script 1.1
# Few comments: 
# - script can submit many many jobs. You need to monitor it. 
# - Currently doesn't take into account the mysql servers load - but the analyses that load the servers. 
# - script consider number of connections to mysql servers. 

#!/usr/bin/python



import sys
import argparse
import logging
import subprocess
import sqlalchemy as db
import time

# from subprocess import Popen, PIPE
from sqlalchemy.sql import select
from sqlalchemy import func, and_
from urllib.parse import urlparse
from sqlalchemy.orm import sessionmaker
from colorama import init
from colorama import Fore, Back, Style
from datetime import datetime

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql
from numpy.f2py.crackfortran import verbose
from test.test_pprint import set2
pymysql.install_as_MySQLdb()



# this is a pipeline object
class Pipeline:
    
    # Constructor
    def __init__(self, name, pipeline_status = 'stop', registry_file = '', mysql_status=0, farm_status=0, pipeline_input_ids_status=0):
        self.name            = name  
        self.pipeline_status = pipeline_status
        self.registry_file   = registry_file
        self.mysql_status    = mysql_servers_status(self.name, 2, verbose=False)
        self.farm_status     = check_farm_jobs() 
        self.pipeline_input_ids_status = check_pipeline_inputIDs_status(self.name, verbose= False)
            
    # Print all info object holds
    def print_data(self, verbose=False):
        if verbose:
            print('# name, {}'.format(self.name))
            print('# pipeline status, {}'.format(colour_pipeline(self.name, self.pipeline_status)))
            print('# mysql_status, {}'.format(self.mysql_status))
            print('# Can I submit more?, {}'.format(self.farm_status))
            print('# pipeline_input_ids_status, {}'.format(self.pipeline_input_ids_status))
    
    # we don't need this one as it is. We used it like: current_status = g.get_set_status()
    def get_set_status(self, new_status = ''):
        if new_status: 
            # print('old_status: ', self.pipeline_status) 
            self.pipeline_status = new_status 
            # print('AND I WILL UPDATE IT TO new_status: ', self.pipeline_status)
        else: 
            # print('just_report_status: ', self.pipeline_status)
            return self.pipeline_status 
    
    # Calculate parameters that will be used by beekeeper. 
    def calculate_bkeeper_params(self, number_of_loops, sleep_div, number_of_workers, verbose = False): 
        # max_of_mysql_connections = sum(self.mysql_status.values()) # TODO: this is not working fine: mysql_status, {'max': 1406}
        max_of_mysql_connections  = mysql_servers_status(self.name, check=2, verbose=False)['max']

        if verbose: 
            print('# script will do the calculations')
            print('max_of_mysql_connections:', max_of_mysql_connections)

        ### Security CHECKS ### 
        # STOP/GO checks
        while max_of_mysql_connections  > 2300:
            print('max_connections are: very high (', max_of_mysql_connections , ')', '. You need to wait.' )
            time.sleep(60)
            max_of_mysql_connections  = sum(mysql_servers_status(self.name, check=2))
        
        while not self.farm_status:
            print('total_number_of running jobs (in farm) are: Very high. You need to wait.' )
            time.sleep(60)
            self.farm_status = check_farm_jobs() 
            # keep sending jobs else sleep.     
        
        run_jobs = check_pipeline_inputIDs_status(self.name, verbose=False).get("RUN", 0)        
        ready_jobs = check_pipeline_inputIDs_status(self.name, verbose=False).get("READY", 0)  
        ready_test = self.pipeline_input_ids_status.get("READY", 0 )

        
        if  ( ( self.pipeline_status == 'boost_run') and ( ready_jobs > 1000 ) ): 
            # this is going to be a real boost for parts of high priority pipelines that makes connections to mysql servers.  
            print('I am in a boost mode here. I will take everything available dude')
            max_loops = 15
            sleep_div = 0.2            
            how_much_space = 2500 - max_of_mysql_connections
            number_of_workers = how_much_space/2 # how_much_space/2 is the high-priority bonus workers
            print('I have: ', how_much_space, ' to use and I will use: ',  number_of_workers)
        elif ( (max_of_mysql_connections < 500 ) and (ready_jobs > 10000 ) ):
            print('+10000_ready_jobs_are_waiting')
            number_of_loops = 15
            sleep_div = 0.2
            number_of_workers = number_of_workers*1.8
        elif ( (max_of_mysql_connections < 1000 ) and ( ready_jobs > 10000 ) ):
            print('Many mysql_connections and +10000_ready_jobs_are_waiting')
            number_of_loops = 8
            sleep_div = 0.3
            number_of_workers = number_of_workers*1.3
        elif ( (max_of_mysql_connections < 1000 ) and ( ready_jobs > 1000 ) ): 
            print('normal')
            number_of_loops = 5
            sleep_div = 0.1 
            number_of_workers = number_of_workers*1.1            
        elif ( (max_of_mysql_connections < 1000 ) and ( ready_jobs < 5 ) 
               and ( ready_jobs > 0 ) ):
            print('not much to run, I will loop twice')
            number_of_loops = 2
        elif ( (max_of_mysql_connections < 1000 ) and ( ready_jobs) ):
            number_of_loops = 10
            print('Many sql connections and very few ready_jobs')
        elif ( ready_jobs == 0 ) :
            if ( run_jobs == 0 ) : 
                print('Nothing_to_run. Looks like we might have finish! ')
                number_of_loops = 1                
        else: 
            print('Problem: Looks like mysql servers are loaded. ')
        
        number_of_workers = int(round(number_of_workers)) 

        # you don't want to submit too many jobs per pipeline
        if number_of_workers > 1100 : 
            number_of_workers = 850
        
        print('n_loops:', number_of_loops, 'sleep:', sleep_div, 'n_workers: ', number_of_workers)
        return number_of_loops, sleep_div, number_of_workers


# Find analyses that are ready to go    
def ready_analyses(db_url, status = 'READY', verbose=False, return_type = 'name'):
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
            # print('ready analysis:',analysis_id, 'logic_name', logic_name)
    else:
        raise ValueError(f'Could not connect to the target database {db_url}.')
    return connections_info


# run commands: I sent all commands and execute them here. 
def block_checkcall(cmd, verbose=False) :
    verbose = 1
    if verbose:
        print('running cmd: ', cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if verbose == 3: 
        print("err:", stderr)
        print("out:", stdout)
    return stdout, stderr 


# check run/pending/fail jobs: 
def check_farm_jobs(task_max=2500, verbose=True):
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


# checking the processlist of mysql servers 
def mysql_servers_status(db_url, check=1, verbose=False):
    url_attributes = urlparse(db_url)
    connections_info = {}

    if check == 1:
        list_cmd = [url_attributes.hostname , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
        cmd_processlist = ''.join(list_cmd)
        stdout_pend, stderr_pend = block_checkcall( cmd_processlist, verbose= False )
        connections_info['total'] = stdout_pend
    elif check == 2:
        # This will calculate max connections across servers
        set1 = ('mysql-ens-genebuild-prod-2','mysql-ens-genebuild-prod-3', 'mysql-ens-genebuild-prod-4')
        set2 = ('mysql-ens-genebuild-prod-5','mysql-ens-genebuild-prod-6', 'mysql-ens-genebuild-prod-7')
        max_connections = int(1) # this value will report the max number of connections. 
        user_time = int(0)
        if url_attributes.hostname == 'mysql-ens-genebuild-prod-4' : 
            for test_host in set1:
                start_time = time.time()
                list_cmd = [test_host, "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
                cmd_processlist = ''.join(list_cmd)
                stdout_pend, stderr_pend = block_checkcall( cmd_processlist, verbose= False )
                stdout_pend = int(stdout_pend)
                if max_connections < stdout_pend: 
                    max_connections = stdout_pend
                
                user_time = time.time() - start_time
                print(f"{test_host} takes : {user_time} sec to respond")
                die_if_server_dont_respond(user_time)    
        elif url_attributes.hostname == 'mysql-ens-genebuild-prod-7' :
            for test_host in set2:
                start_time = time.time()
                list_cmd = [test_host , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | " ,
                     " wc -l" ]
                cmd_processlist = ''.join(list_cmd)
                stdout_pend, stderr_pend = block_checkcall( cmd_processlist, verbose= False)
                stdout_pend = int(stdout_pend)
                if max_connections < int(stdout_pend): 
                    max_connections = stdout_pend 
                user_time = time.time() - start_time
                print(f"{test_host} takes : {user_time} ") 
                die_if_server_dont_respond(user_time)
        else :
            print('problem')

        connections_info['max'] = max_connections

    else: 
        list_cmd = [ url_attributes.hostname , "-ensadmin " , " -e ", " 'SHOW PROCESSLIST' ", " | ", 
                    " awk ", " '{print $5}' "  ]
        cmd_processlist = ''.join(list_cmd)
        stdout_pend, stderr_pend = block_checkcall( cmd_processlist, verbose= False )

        for ser_status in stdout_pend.split(): 
            if (ser_status in connections_info): 
                connections_info[ser_status] += 1 
            else: 
                connections_info[ser_status] = 1 

    if verbose: 
        print("Verbose output of mysql_servers_status def")
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
    if (pipeline_stage in ['run', 'boost_run'] ): 
        print(Fore.GREEN + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'stop'):
        print(Fore.MAGENTA + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'FAILED'):    
        print(Fore.RED + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'ready'):
        print(Fore.YELLOW + pipeline_name , pipeline_stage )
    elif (pipeline_stage == 'done'):
        print(Fore.WHITE + pipeline_name , pipeline_stage )
    else: 
        print('are you sure about the status?') 

        
# reset boost priority 
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
        analyses_mysql_free = ['download_long_read_fastq', 'minimap2', 'minimap2_himem', 'download_RNASeq_fastqs','bwa',\
                     'bwa2bam','merged_tissue_file','bam2bigwig', 'create_analyses_type_job']
        return analyses_mysql_free
    elif list_type == 'mysql_Load': 
        # analyses that load mysql servers (experience speaking)
        analyses_mysql_load = ['wga_project_transcripts','layer_annotation','run_utr_addition','genebuilder', 'rebatch_repeatmasker', 
                     'run_repeatmasker_small_batch', 'remove_redundant_genblast_genes']
        return analyses_mysql_load
    elif list_type == 'quickly_done':
        # analyses that are quickly done and can be loop more times when in high priority
        # https://docs.google.com/spreadsheets/d/1RgPW9q05YGENajVs9XsEr_XN7qZCM7rGt-_d0z3iZ1Y/edit?usp=sharing
        analyses_that_run_quickly = ['core_assembly_name_update', 'otherfeatures_assembly_name_update', 'pairaligner_stats',\
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
                                      'change_biotype_for_weak_cds']
        return analyses_that_run_quickly
    else:
        print('problem')    


# This function will check for long running jobs, usually there is a problem when a job runs for ever.  
def check_for_long_running_jobs():
    return


# connect to registry db and update status of the pipeline

# display basic table: 
def create_html(data):
    # https://www.ebi.ac.uk/~kbillis/test/display/
    # current date and time
    dateTimeObj = datetime.now()

    new_data = {}
    for i in data:
        
        url_attributes = urlparse(i)
        # url_attributes = urlparse(i)
        info_to_print = url_attributes.hostname + url_attributes.path

        i_new = i
        # i_new = i_new.replace('ensembl','xxxxxxx')
        new_data[info_to_print] = data[i].get_set_status()
        # print("old:",i,"new:",i_new, "data_old:", data[i],"data_new:", new_data[i_new])
    data = new_data    

    html = '<table><tr><th>' 
    for row in data:
        html += '<tr><td>' + str(row) +  '</td><td>'  + str(data[row]) + '</td><td>' + str(dateTimeObj) + '</td></tr>' 
    html += '</table>'

    file_html = '/homes/kbillis/public_html/test/display/index.html'
    with open(file_html, 'w') as newHTMLFile:
        newHTMLFile.write(html)
        newHTMLFile.close()


def die_if_server_dont_respond(respond_time):
    if (respond_time >10): 
        sys.exit('MAJOR ISSUE WITH MYSQL SERVER ! emergency stop')
    elif (respond_time > 3): 
        print('Respond time of mysql servers is high, I will sleep a bit. ')
        time.sleep(60)
    elif (respond_time > 1): 
        print('Respond time of mysql servers is not good ')
    else: 
        print('Respond time of mysql servers is OK')



def main():
    """
    main function
    """
    
    # colouring names. 
    init(autoreset=True)

    # get arguments from command line: 
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--input", help="Your input file. Example: /hps/nobackup2/production/ensembl/kbillis/production/genebuilds/2020_11_python_test/running_pipelines.txt")
    # parser.add_argument("-r", "--reg_conf", help="Your registry/Databases.pm file if you are ")
    parser.add_argument("-a", "--act", help="What I have to do?")
    parser.add_argument("-l", "--number_of_loops", type=int, help="How many loops?", default=10 )
    parser.add_argument("-w", "--number_of_total_workers", type=int, help="How many total_workers?", default=250 )
    parser.add_argument("-s", "--skip_checks", help="Skip tests?", default=False, action='store_true')
    parser.add_argument("-b", "--skip_boost", help="Skip priority?", default=False, action='store_true')
    parser.add_argument("-v", "--verbose",dest='verbose',action='store_true', help="Verbose mode.")
    parser.add_argument("-p", "--analyses_pattern_str", help="Any analyses pattern to run?")    
    
    options = parser.parse_args()
    verbose = options.verbose
    debug_parameter = '  '
    if verbose:
        print("Verbose mode on")
        debug_parameter = ' -debug 1 '
    else:
        print("Verbose mode off")



    # standard parameters: 
    beekeeper_locations = '/nfs/production/panda/ensembl/kbillis/enscode_2021_03/enscode/ensembl-hive/scripts/beekeeper.pl '
    # debug_parameter = ' -debug 1 '
    debug_parameter = ' '
    how_many_loops = options.number_of_loops
    priority_limit = int(how_many_loops*0.1)
    skip_checks = options.skip_checks
    skip_boost = options.skip_boost
    analyses_pattern = options.analyses_pattern_str
    # registry_file = options.reg_conf

    if not options.act == 'loop': 
        skip_boost = 1 
        how_many_loops = 1
    
    ### Start running the script: ### 
    # read a dictionary with pipelines to run: 
    d = {}
    try: 
        with open(options.input) as f:
            # data = f.readlines()
            for line in f:
                pipeline_status = 'empty'
                if line.startswith('#'):
                    print('WARN: You comment out this pipeline: ',line)
                else:
                    line_info = line.rstrip()
                    pl_name,pl_registry = line_info.split(' ')
                    if not skip_boost : 
                        boost_option = input('Do you want to boost pipeline ' + pl_name + '?(y/n)')
                        if boost_option.startswith('y'): 
                            pipeline_status = 'boost_run'
                        else: 
                            pipeline_status = 'run'
                    else:
                         pipeline_status = 'run'
    
                    print('#LOADING# pipeline to start: ', pl_name, ' the status is: ', pipeline_status, 'and the Database file is:', pl_registry)
                    d[pl_name] = Pipeline(pl_name, pipeline_status, pl_registry)
    except OSError as ex: 
        print('Error: File is not available') 

    
    # this could be a while: 
    for ii in range(how_many_loops): 
        # if ii == priority_limit:
        #    d = reset_priorities(d)
        print('\n###########################Start of loop###############################')
        
        # check how many workers each pipeline should run. 
        bkeeper_total_workers = options.number_of_total_workers
        worker_equal_dis = int(float(bkeeper_total_workers/len(d) ) ) 
        
        # this is just informative and useful if all pipelines die: 
        amount_of_failed_pipelines = 0 
        for pipeline_name in d: 
            g = d[pipeline_name]
            current_status = g.pipeline_status
            if (current_status == 'FAILED'):
                amount_of_failed_pipelines = amount_of_failed_pipelines +1
          

        if amount_of_failed_pipelines == len(d): 
            sys.exit('EXIT MESSAGE: ALL RUNNING PIPELINES FAILED ! ')
        elif amount_of_failed_pipelines > 0 : 
            print(f"WARNING: You have {amount_of_failed_pipelines} failed pipelines.")
        else: 
            print(f"INFO: You have 0 FAILED PIPELINES and {len(d)} pipelines to run.")

            
        ## loop those pipelines 
        for pipeline_name in d:
            # Get an instance of the Pipeline class
            print('#########PIPELINE ', pipeline_name, 'START #########')
            g = d[pipeline_name]
            
            bkeeper_max_loops = 5
            bkeeper_sleep = 0.3

            g.pipeline_input_ids_status = check_pipeline_inputIDs_status(pipeline_name, verbose= False)
            current_status    = g.pipeline_status
            g.mysql_status    = mysql_servers_status(pipeline_name, 2, verbose=False)
            g.farm_status     = check_farm_jobs()

            g.print_data(verbose=True)
            
            if not options.act == 'loop':
                mylist = ['perl ', beekeeper_locations, ' ', debug_parameter , ' -' , options.act , ' -url ' , pipeline_name ]
                cmd_tmp = ''.join(mylist)
                stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose)
                # Construct an instance of the Pipeline class   
                # g.print_data(verbose=True)   # Call an instance method; prints 
            else: 
                if current_status == 'FAILED' :
                    print(pipeline_name, ' pipeline has FAILED issues:')
                    colour_pipeline(pipeline_name, current_status) 
                    
                    # Are there still failed_jobs? : 
                    failed_jobs = check_pipeline_inputIDs_status(pipeline_name, verbose=False).get("FAILED", 0)
                    if failed_jobs == 0:
                        print('No failed jobs for this pipeline.')
                        g.pipeline_status = 'run'
                    else:
                        print('Failed jobs for this pipeline. Number of failed_jobs: ', failed_jobs)
                        continue

                    
                
                elif current_status == 'config' :
                    # method to create config
                    print(pipeline_name, 'I will create config')
                    g.pipeline_status = 'run'
                    continue
                
                # SET beekeeper parameters based on analysis you are running/read_to_run
                analyses_ready = ready_analyses(pipeline_name)
                print('# Analyses ready (status) to run/running: ', '[%s]' % ', '.join(map(str, analyses_ready)))
                analyses_mysql_free = mysql_load_analyses(list_type='mysql_Free')
                is_mysql_Free = False
                is_mysql_Free = all(x in analyses_mysql_free for x in analyses_ready)

                analyses_mysql_heavy = mysql_load_analyses(list_type='mysql_Load')
                is_mysql_Heavy = False
                is_mysql_Heavy = any(x in analyses_mysql_heavy for x in analyses_ready)


                # include Databases.pm as parameter or not. 
                registry_filename = g.registry_file
                analyses_ready_analysis_ids = analyses_ready.values()
                analysis_id_to_run = 0
                if analyses_ready_analysis_ids: 
                    analysis_id_to_run = max(analyses_ready_analysis_ids)
                
                registry_param = ''
                print("DEBUG:: analysis_to_run::", analysis_id_to_run)
                if ( (registry_filename ) and (analysis_id_to_run > 70 ) ) :   # if you are above 40, I assume you created Databases.pm
                    try:
                        my_file = open(registry_filename)
                    except IOError as err:
                        print(err.errno)
                        print(err.strerror) 
                       
                    registry_param = ' -reg_conf ' + str(registry_filename)
                else:
                    print("No registry/Databases.pm provided - this might cause issues in lastZ part")

                analyses_pattern_param = ''
                if (analyses_pattern): 
                    analyses_pattern_param = ' -analyses_pattern ' + ' "' + analyses_pattern + '" '
                
                if skip_checks: 
                    bkeeper_total_workers = worker_equal_dis
                elif is_mysql_Heavy: 
                    bkeeper_total_workers = int(worker_equal_dis/3) # this is a punishment, but I can't let it run too many jobs. 
                    bkeeper_max_loops = 8
                    bkeeper_sleep = 0.2
                    print('# You are running a mysql HEAVY analysis.  workers are: ', bkeeper_total_workers, ', bkeeper_max_loops:'
                          , bkeeper_max_loops , ' and bkeeper_sleep: ', bkeeper_sleep)                
                elif is_mysql_Free: 
                    bkeeper_total_workers = 500
                    bkeeper_max_loops = 8
                    bkeeper_sleep = 0.2
                    print('# You are running a mysql FREE analysis.  workers are: ', bkeeper_total_workers, ', bkeeper_max_loops:'
                          , bkeeper_max_loops , ' and bkeeper_sleep: ', bkeeper_sleep)
                else:
                    # Here is the semi-dynamic part
                    # Construct an instance of the Pipeline class
                    print('# You are running a mysql NORMAL analysis. workers are: ', bkeeper_total_workers, ', bkeeper_max_loops:'
                          , bkeeper_max_loops , ' and bkeeper_sleep: ', bkeeper_sleep)

                    # g.print_data(verbose=True)  
                    # will return number of workers, amount of loops and sleep value.
                    bkeeper_max_loops, bkeeper_sleep, bkeeper_total_workers = g.calculate_bkeeper_params(bkeeper_max_loops, bkeeper_sleep, worker_equal_dis)

                # build beekeeper to run: 
                mylist = ['perl ' , beekeeper_locations , ' ' , debug_parameter , ' -' , options.act , ' -url ' , pipeline_name ,\
                          ' -max_loops ' , str(bkeeper_max_loops) , ' -total_running_workers_max ' , str(bkeeper_total_workers) ,\
                          ' -sleep ' , str(bkeeper_sleep) , registry_param, analyses_pattern_param]
                cmd_tmp = ''.join(mylist)
                stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose=False)

                # this is a high-priority and it will run until jobs fail or if there are analyses with single jobs that doesn't require usually more than 2 sec
                if (current_status == 'boost_run'):
                    ready_jobs = check_pipeline_inputIDs_status(pipeline_name, verbose=False).get("READY", 0)
                    print('Extra loops: I like this pipeline (boost run), ', pipeline_name, 'and I will run again and again until ready jobs are low.')
                    while ready_jobs > 5000: 
                        # give it a bit more of life until more input_id are done.
                        bkeeper_max_loops, bkeeper_sleep, bkeeper_total_workers = g.calculate_bkeeper_params(bkeeper_max_loops, bkeeper_sleep, worker_equal_dis)
                        mylist = ['perl ' , beekeeper_locations , ' ' , debug_parameter , ' -' , options.act , ' -url ' , pipeline_name ,\
                                  ' -max_loops ' , str(bkeeper_max_loops) , ' -total_running_workers_max ' , str(bkeeper_total_workers) ,\
                                  ' -sleep ' , str(bkeeper_sleep) , registry_param, analyses_pattern_param]
                        cmd_tmp = ''.join(mylist)
                        stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose) 
                        ready_jobs = check_pipeline_inputIDs_status(pipeline_name, verbose=False).get("READY", 0)
                    
                    analysis_list = mysql_load_analyses(list_type='quickly_done')
                    while analyses_ready in analysis_list:
                        analyses_ready = ready_analyses(pipeline_name)
                        stdout_run, stderr_run = block_checkcall(cmd_tmp, verbose) 
                
                
                # if there are long running jobs this update pipeline to something else: 
                # major check if there is an analysis that takes too much time (for example, if one analysis runs for more than 10 hours,
                # if failed jobs this will be change the status 
                check_failed_j = 1              
                failed_jobs = check_pipeline_inputIDs_status(pipeline_name, verbose=False).get("FAILED", 0)
                if failed_jobs == 0:
                    print('No failed jobs for this pipeline.')
                else:
                    print('Failed jobs for this pipeline. Number of failed_jobs: ', failed_jobs)
                    if check_failed_j > 0:
                        g.pipeline_status = 'FAILED'
                        
            print('#########PIPELINE ', pipeline_name, ' END #########')



        dateTimeObj = datetime.now()
        print('last update: ', dateTimeObj)  

        # TODO : check for many things 
        
        # create html that will report things. 
        create_html(d)
        # call the class 
        # pipeline, mysql_check, check_jobs 
        print('\n###########################End of', ii ,' loop###############################')
        
    print('Game over!')




if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print()
        print("Interrupted with CTRL-C, exiting...")
        sys.exit()

    
