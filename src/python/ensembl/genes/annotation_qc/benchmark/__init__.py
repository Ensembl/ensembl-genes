# Copyright 2026 EMBL-EBI
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Reproducible annotation-comparison experiments.

An experiment configuration names references, predictor runs and evaluation
policies. This package validates and prepares inputs, runs (or imports)
``pairwise-compare`` results with content-based caching, identifies the code that
produced each result (registry.py), and builds the SQLite dataset that the
dashboard reads by indexing the comparator's contract outputs verbatim (pre-contract
results only through the explicit legacy importer). It computes no metric:
comparison semantics live in ``metrics/pairwise`` and are not changed here.
"""
