#!bash
(
# Exit on error
set -o errexit

# Set up the directories for data
mkdir database
cd database
mkdir differentiability group_action linear_subspaces zeta zeta_functions

cd differentiability
mkdir singular smooth

cd ../group_action
mkdir orbit_representatives stabilizers_info
mkdir orbit_representatives/orbit_representative_in_V

cd ../linear_subspaces
mkdir lines_through_cubics planes_through_cubics smooth_lines_singular_cubics

# Download and unpack the cubic orbit.
cd ..
cd group_action/orbit_representatives/orbit_representative_in_V

for i in $(seq 1 85);
do
    curl https://arxiv.org/src/2306.09908v1/anc/ancillary/database/group_action/orbit_representatives/orbit_representative_in_V/orbitreps-$i.data \
         -o orbitreps-$i.data
done

)
