
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

import time
from builtins import input


class hello_world(eHive.BaseRunnable):
    """Runnable to Hello"""

    def param_defaults(self):
        return {
            'take_time' : 0
        }


    def run(self):
        a_multiplier = self.param_required('hello')
        print("Hi I am Python. " + str(a_multiplier))
        time.sleep( self.param('take_time') )
        #input_check = input("Do you want me to check your db?" )
        #if (input_check) : 
        #    print("you said: " + str(input_check))


    def write_output(self):
        self.dataflow({},1)
    
