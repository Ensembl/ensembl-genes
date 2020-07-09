
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
import mysql.connector
from mysql.connector import errorcode


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
        alist=["mysql-ens-genebuild-prod-6","4532", "ensro", "kbillis_sarcophilus_harrisii_core_100"]
        connect_to_db(alist)


    def write_output(self):
        self.dataflow({},1)
    
    

def connect_to_db(db_connection_param):
        
    print("And all the rest... %s" %(list(db_connection_param)))
    # for item in db_connection_param:
    for position in range(len(db_connection_param)):
        print(db_connection_param[position])
        print(position)
    
    server_read,port_read,user_read,database_read = db_connection_param

    print("server::" + server_read)

    status = "good"
    try:
        cnx = mysql.connector.connect(user=user_read, database=database_read, port=port_read, password="", host=server_read)
        # cnx = mysql.connector.connect(user='ensro', database='', port=4532, password='', host='mysql-ens-genebuild-prod-6')
        # cnx = mysql.connector.connect(user='ensro', port=4532, password='', host='mysql-ens-genebuild-prod-6', buffered=True)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
        else:
            print(err)
    else:
        print("all ok")
        
    # what to do: check all the dbs
    todo = "query" 
    if todo == "show" :
        cursor = cnx.cursor()
        databases = ("show databases")
        cursor.execute(databases)
        for (databases) in cursor: 
            print(databases[0])
    else: 
        cursor = cnx.cursor()
        query = ("SELECT * FROM meta")
        cursor.execute(query)
        remaining_rows = cursor.fetchall()
        # print(cursor)
        for (i) in range(len(remaining_rows)):
            print(i," is ", remaining_rows[i] )
        
    cursor.close()
    cnx.close()
    return status
