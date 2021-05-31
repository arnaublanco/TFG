### fMRIPrep Processing Pipeline ###

bids_root_dir=$HOME/Desktop/ds002715 # Root folder
nthreads=4 # Maximum number of threads
mem=20 # Upper bound RAM memory limit
mem=`echo "${mem//[!0-9]/}"` # Remove gb suffix at the end (for fmriprep-docker command)
mem_mb=`echo $(((mem*1000)))` # Gb -> Mb

# Use FreeSurfer license
export FS_LICENSE=$HOME/Documents/license.txt

# For loop that iterates through all the subjects
for i in {01..08}
do

subj=$i

# Run fmriprep #

# Input: $bids_root_dir
# Output: $bids_root_dir/derivatives

fmriprep-docker $bids_root_dir $bids_root_dir/derivatives \
    participant \
    --participant-label $subj \
    --skip-bids-validation \
    --md-only-boilerplate \
    --fs-license-file $HOME/Documents/license.txt \
    --output-spaces MNI152NLin2009cAsym fsnative \
    --nthreads $nthreads \
    --stop-on-first-crash \
    --mem_mb $mem_mb \
    -w $HOME

done
