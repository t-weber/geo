#!/bin/bash
#
# get external dependencies
# @author Tobias Weber (orcid: 0000-0002-7230-1932)
# @date 1-jun-2021
# @license see 'LICENSE' file
#

mkdir -p ext
pushd ext
wget https://raw.githubusercontent.com/boostorg/polygon/develop/example/voronoi_visual_utils.hpp
popd
