
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

import eHive

import datetime
# import mysql.connector

import time
#from builtins import input


class connect_to_mysql_db(eHive.BaseRunnable):
    """Runnable to connect to db"""

    def param_defaults(self):
        return {
            'take_time' : 0
        }


    def run(self):
        # a_multiplier = self.param_required('hello')
        print("Hi I am Python and connecto to db")
        connect_to_db("mysql-ens-genebuild-prod-6","4532", "ensro", "kbillis_sarcophilus_harrisii_core_100")


    def write_output(self):
        self.dataflow({},1)
    
    

def connect_to_db(*db_connection_param):
        
    status = "OK"
    print("And all the rest... %s" %(list(db_connection_param)))

        # status = "good"
        # cnx = mysql.connector.connect(user='ensro', database='')
        # cursor = cnx.cursor()
        
        # query = ("SELECT first_name, last_name, hire_date FROM employees "
        # "WHERE hire_date BETWEEN %s AND %s")
        
        #hire_start = datetime.date(1999, 1, 1)
        #hire_end = datetime.date(1999, 12, 31)
        
        #cursor.execute(query, (hire_start, hire_end))
        
        #for (first_name, last_name, hire_date) in cursor:
        #    print("{}, {} was hired on {:%d %b %Y}".format(
        #        last_name, first_name, hire_date))
        
        #cursor.close()
        #cnx.close()
    return status
