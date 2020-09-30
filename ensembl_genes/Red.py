# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import filecmp
import errno
import subprocess

import sqlalchemy as db
import sqlalchemy_utils as db_utils

# sqlalchemy requires MySQLdb but MySQLdb doesn't support Python 3.x
# pymysql can be imported and used instead
import pymysql
pymysql.install_as_MySQLdb()

import eHive

class Red(eHive.BaseRunnable):
    """Runnable that runs Red to find repeats and store them in the target database"""

    def param_defaults(self):

        return {
            'logic_name' : 'red',
            'target_db' : {'dbname' : '', \
                           'driver' : '', \
                           'host' : '', \
                           'pass' : '', \
                           'port' : '', \
                           'user' : ''},
            'red_path' : '',
        }

    def fetch_input(self):

        genome_file = self.param_required('genome_file')
        gnm = self.param('gnm',self.param_required('genome_file_tmpdir'))
        msk = self.param_required('msk')
        rpt = self.param_required('rpt')
        red_path = self.param_required('red_path')
        target_db = self.param_required('target_db')
        self.param('target_db_url',target_db['driver']+'://'+ \
                                   target_db['user']+':'+ \
                                   target_db['pass']+'@'+ \
                                   target_db['host']+':'+ \
                                   str(target_db['port'])+'/'+ \
                                   target_db['dbname']+'?'+ \
                                   'local_infile=1')

        # make sure the genome_file tmpdir exists and copy the genome file into it
        # in this way we make sure that the only .fa file to be processed is the one we want
        try:
            os.makedirs(gnm,exist_ok=True)
            os.chmod(gnm,0o777)
        except:
            print('Could not create '+gnm+' directory in "genome_file_tmpdir".')
            raise PermissionError(errno.EPERM,os.strerror(errno.EPERM),gnm)

        new_genome_file = gnm+'/'+os.path.basename(genome_file)
        try:
            os.popen('cp '+genome_file+' '+new_genome_file).read() # without read() Python keeps running before the file is written
        except:
            print('Could not copy file '+genome_file+' into directory '+gnm)
            raise PermissionError(errno.EPERM,os.strerror(errno.EPERM),new_genome_file)

        # check that the file genome_file exists, it ends with .fa
        # and there is not any other .fa file within the same directory
        if sum(map(lambda x: x.endswith('.fa'),os.listdir(gnm))) > 1:
            raise ValueError('The Red program requires that the directory containing the .fa file does not contain any other .fa file. The directory '+gnm+' in "gnm"  contains more than one .fa file.')
        elif not(os.path.isfile(genome_file)):
            raise ValueError('The file '+genome_file+' in "genome_file" does not end with .fa as required by the Red program.')

        # check that the file was copied successfully
        if not(filecmp.cmp(genome_file,new_genome_file,shallow=False)):
            raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),new_genome_file)

        genome_file = self.param(genome_file,new_genome_file)

        # check that the Red binary exists
        if not(os.path.isfile(red_path)):
            raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),red_path)

        # check that the target database exists
        if not(db_utils.database_exists(self.param('target_db_url'))):
            raise ValueError('Could not connect to the target database '+target_db+' in "target_db".')

        # make sure that the output directories exist
        try:
            os.makedirs(msk,exist_ok=True)
        except:
            print('Could not create '+msk+' directory in "msk".')
            raise PermissionError(errno.EPERM,os.strerror(errno.EPERM),msk)

        try:
            os.makedirs(rpt,exist_ok=True)
        except:
            print('Could not create '+rpt+' directory in "rpt".')
            raise PermissionError(errno.EPERM,os.strerror(errno.EPERM),rpt)

        # check that the output directories are empty
        if os.listdir(msk):
            raise ValueError('The msk output directory '+msk+' is not empty.')

        if os.listdir(rpt):
            raise ValueError('The rpt output directory '+rpt+' is not empty.')


    def run(self):

        # output format: 1 (chrName:start-end) or 2 (chrName start end)
        # Note that chrName includes the '>' character
        cmd = self.param('red_path')+ \
              ' -frm 2'+ \
              ' -gnm '+self.param('gnm')+ \
              ' -msk '+self.param('msk')+ \
              ' -rpt '+self.param('rpt')

        try:
            response = subprocess.check_call(cmd.split())
        except subprocess.CalledProcessError as err:
            print("Could not run Red. Return code "+err.returncode)


    def write_output(self):

        engine = db.create_engine(self.param('target_db_url'))
        connection = engine.connect()
        metadata = db.MetaData()

        analysis_table = db.Table('analysis',metadata,autoload=True,autoload_with=engine)
        meta_table = db.Table('meta',metadata,autoload=True,autoload_with=engine)
        repeat_consensus_table = db.Table('repeat_consensus',metadata,autoload=True,autoload_with=engine)
        repeat_feature_table = db.Table('repeat_feature',metadata,autoload=True,autoload_with=engine)

        # insert Red analysis
        analysis_insert = analysis_table.insert().values({'created':db.sql.func.now(), \
                                                         'logic_name':self.param('logic_name'), \
                                                         'program':'Red', \
                                                         'program_version':'05/22/2015', \
                                                         'program_file':self.param('red_path')})
        connection.execute(analysis_insert)

        # fetch the inserted analysis_id
        analysis_query = db.select([analysis_table.columns.analysis_id])
        analysis_query = analysis_query.where(analysis_table.columns.logic_name == 'red')
        analysis_results = connection.execute(analysis_query).fetchall()
        analysis_id = analysis_results[0][0]

        # insert repeat analysis meta keys
        meta_insert = meta_table.insert().values({'species_id':1, \
                                                  'meta_key':'repeat.analysis', \
                                                  'meta_value':'red'})
        connection.execute(meta_insert)

        # insert dummy repeat consensus
        repeat_consensus_insert = repeat_consensus_table.insert().values({'repeat_name':'Red', \
                                                                          'repeat_class':'Red', \
                                                                          'repeat_type':'Red', \
                                                                          'repeat_consensus':'N'})
        connection.execute(repeat_consensus_insert)

        # fetch the inserted repeat_consensus_id
        repeat_consensus_query = db.select([repeat_consensus_table.columns.repeat_consensus_id])
        repeat_consensus_query = repeat_consensus_query.where(repeat_consensus_table.columns.repeat_name == 'Red')
        repeat_consensus_results = connection.execute(repeat_consensus_query).fetchall()
        repeat_consensus_id = repeat_consensus_results[0][0]

        # parse the repeats file and make a tsv file ready to load into the repeat_feature table
        repeats_file = self.parse_repeats(self.param('rpt'),repeat_consensus_id,analysis_id)

        # insert repeat features
        repeat_feature_query = "LOAD DATA LOCAL INFILE '"+repeats_file+"' INTO TABLE repeat_feature FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n' (seq_region_id,seq_region_start,seq_region_end,repeat_start,repeat_end,repeat_consensus_id,analysis_id)"
        repeat_feature_result = connection.execute(repeat_feature_query)


    def parse_repeats(self,rpt,repeat_consensus_id,analysis_id):

        rpt_files = os.listdir(rpt)
        if not rpt_files:
            raise FileNotFoundError(errno.ENOENT,os.strerror(errno.ENOENT),"The repeats output directory "+rpt+" is empty.")
        elif len(rpt_files) > 1:
            raise ValueError('The rpt output directory '+rpt+' contains more than 1 file.')
        elif not rpt_files[0].endswith('.rpt'):
            raise ValueError('The file '+rpt_files[0]+' in the rpt output directory '+rpt+' does not end with .rpt as required.')
        else: # 1 file in rpt dir and ends with .rpt as expected
              # Red's rpt output file contains ">" which needs to be removed from each line
              # and we need to replace the seq region name with seq region id and add some extra columns so it can be loaded directly
            engine = db.create_engine(self.param('target_db_url'))
            connection = engine.connect()
            metadata = db.MetaData()

            sr = db.Table('seq_region',metadata,autoload=True,autoload_with=engine)
            sra = db.Table('seq_region_attrib',metadata,autoload=True,autoload_with=engine)
            at = db.Table('attrib_type',metadata,autoload=True,autoload_with=engine)

            query = db.select([sr.columns.seq_region_id,sr.columns.name])
            query = query.select_from(sr.join(sra, \
                                              sr.columns.seq_region_id == sra.columns.seq_region_id). \
                                         join(at,
                                              at.columns.attrib_type_id == sra.columns.attrib_type_id))
            query = query.where(at.columns.code == 'toplevel')

            results = connection.execute(query).fetchall()

            seq_region = {}
            for seq_region_id,name in results:
                seq_region[name] = seq_region_id

            rpt_file = rpt+'/'+rpt_files[0]
            fixed_rpt_file = rpt_file+'.fixed'

            with open(rpt_file,'r') as f_in,open(fixed_rpt_file,'w') as f_out:
                for line in f_in:
                    columns = line.split()

                    if columns[0][0] == ">":
                        name = columns[0][1:] # remove first character '>'
                    else:
                        name = columns[0]

                    seq_region_start = int(columns[1])+1 # Red's start is zero-based
                    seq_region_end = int(columns[2])-1   # Red's end is exclusive
                    seq_region_id = seq_region[name]
                    # seq_region_id seq_region_start seq_region_end repeat_start repeat_end repeat_consensus_id analysis_id
                    print(str(seq_region_id)+"\t"+ \
                          str(seq_region_start)+"\t"+ \
                          str(seq_region_end)+"\t"+ \
                          "1\t"+ \
                          str(seq_region_end-seq_region_start+1)+"\t"+ \
                          str(repeat_consensus_id)+"\t"+ \
                          str(analysis_id),
                          file=f_out)
                    #f_out.write(line[1:])

        return fixed_rpt_file
