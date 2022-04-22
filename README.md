# Trapping of Internal Tides in a Channel Model

Two tutorials of the Utrecht University course NS-MO447M Waves in Geophysical Fluids will be devoted to recovering and extending the result published in the paper

> Drijfhout, S., & Maas, L. R. M. (2007). Impact of Channel Geometry and Rotation on the Trapping of Internal Tides. *Journal of Physical Oceanography*, *37*(11), 2740–2763. https://doi.org/10.1175/2007JPO3586.1

[which you can download here](https://webspace.science.uu.nl/~maas0131/files/drijfhoutmaas07jpo%5bsmallpdf.com%5d.pdf).

You will do so using the simple ocean circulation model named “MICOM” (Miami Isopycnic Coordinate Ocean Model Output), written in Fortran, which is the main computer language for nearly all ocean, atmosphere, and climate models.

This repository holds the model code and files that you can use for analysis. This README also contains the instructions for a small report that forms the fifth assignment of the course



## Model description

The model is configured as a channel with an open boundary at one, oceanward side, and a continental slope at the opposite side. The length of the channel is 1200km, the width 191.25 km. There are 43 layers of 100-m depth, so that total depth in the middle of the channel is 4300 m. 

A barotropic tidal wave, having a typical period of 12 h, enters the channel from the open boundary. Over the continental slope its mainly horizontally-moving barotropic tidal current aquires a vertical component, as the flow is forced to follow the bottom. This vertical velocity displaces isopycnals up and downward out-of-equilibrium. The gravitational restoring force subsequently generates outward-propagating internal tides. 

The paper discusses 4 experiments. In two experiments the earth’s rotation is neglected and the Coriolis frequency is set to 0. In the other two experiments we use an <img src="https://render.githubusercontent.com/render/math?math={f}"> plane with <img src="https://render.githubusercontent.com/render/math?math={f = 10^{-4}\ \rm{s}^{-1}}">. Both cases are run with two channel configurations: in one case the channel has a flat bottom and vertical walls in the cross-channel direction; in the other case the channel has a parabolic cross-channel bottom profile with maximum depth of 4300 m and sloping sidewalls. When <img src="https://render.githubusercontent.com/render/math?math={f=0}"> the buoyancy frequency is chosen as <img src="https://render.githubusercontent.com/render/math?math={N = 3.05 \times 10^{-3}\ \rm{s}^{-1}}">. When <img src="https://render.githubusercontent.com/render/math?math={f = 10^{-4}\ \rm{s}^{-1}}">,  <img src="https://render.githubusercontent.com/render/math?math={N= 2.2 \times 10^{-3}\ \rm{s}^{-1}}">. 

The model files in this repository are configured for the case <img src="https://render.githubusercontent.com/render/math?math={f = 10^{-4}\ \rm{s}^{-1}}">, <img src="https://render.githubusercontent.com/render/math?math={N = 2.2 \times 10^{-3}\ \rm{s}^{-1}}">. 



##  Downloading and building the model

You are used to working with Python, which is an interpreted programming language. This means that your instructions are not immediately understood by your computer, but are interpreted and translated to machine code on-the-fly as you execute a script. In contrast, Fortran is a compiled programming language. Before you run a program, your entire code needs to be translated to machine code first. This extra step allows your code to run incredibly efficiently whenever it is executed.

To run the model on your machine, you will first need to compile it. This is also called *building* the model. The model should be able to build on MacOS and Linux with the Fortran compiler *gfortran* (part of the *GCC* compiler collection) installed. We have also tested the model on Gemini, the Science department’s computer cluster. If you are using MacOS, make sure to have the Command Line Tools installed: open the Terminal and type `xcode-select --install`. If you are unable to build the model on your own computer, feel free to use the Gemini cluster. In that case, make sure to read the specific instructions below.

<details>
  <summary>⚠️ Instructions for using the Gemini cluster</summary>

  ### Logging in
  1. Open a Terminal.
  2. Connect to the Gemini cluster by typing `ssh 1234567@gemini.science.uu.nl` using your Solis-ID in place of 1234567.
  3. Type your Solis-ID password.
  4. You're in! Your home directory is `/nethome/1234567`. It has a quotum of 2GB. For temporarily storing large amounts of data, create a scratch folder on the scratch disk: `mkdir /scratch/1234567`. Gemini consists of several computing nodes. Each node has its own `scratch` disk.  This means you may have to copy files between nodes. The default (login) node is `science-bs35`. The other node for regular, scheduled jobs is `science-bs37`. You can switch to this node by typing `ssh science-bs37`, and exit it again by typing `exit`. 

  ### Executing commands
Now you can execute shell commands as usual. Note that Gemini is a cluster, which means that you are sharing resources with other users. When executing small jobs (e.g. copying files, running small scripts, building the model), you can do so as usual. **However, when running larger jubs such as running the model, you should make use of the queueing system.** This way, you're getting adequate computational resources for your _job_, and by using a queue, you won't be hogging resources from other users. 

  ### Submitting a job to the queue
  1. To tell the queueing system what your _job_ is comprised of, you should first create a *job script*, e.g. `my_job.sh`. The `sh` extension indicates that this is a shell script. An example job script looks as follows
  ```bash
  #/bin/bash

  # SGE: Set the job name
  #$ -N wgf_model
  # SGE: this flag exports all active environment variables to the job
  #$ -V
  # SGE: time limit and queue. You probably don't need to change this
  #$ -l h_rt=2:00:00 
  #$ -q all.q
  # SGE: your Email here, for job notification
  #$ -M my.student.email@uu.nl
  # SGE: when do you want to be notified (b : begin, e : end, s : error)?
  #$ -m e 
  #$ -m s
  # SGE: ouput in the current working dir
  #$ -wd /scratch/1234567/wgf_model/

  # Navigate to the right directory and run the model
  cd /scratch/1234567/wgf_model/ # Navigate to where the model is stored
  ./micom1.x # This line starts the model execution
  
  ```
  2. Submit the job using `qsub /path/to/my_job.sh`.
  3. You can inspect the job status using `qstat.`
  4. If you need to delete the job, check the id using `qstat` and use `qdel 123` with 123 being the job id.

  ### Running Jupyterlab on the cluster
  You can use Jupyter Lab on the cluster. This allows you to easily analyze the model output. 
  1. To do so, you must first load _Conda_: `module load miniconda/3`. Initialize Conda by typing `conda init bash`. You may need to open another bash-shell: type `bash`. You can tell that Conda is loaded when `(base)` is being shown in front of the interpreter.
  2. Start Jupyter: `jupyter lab --no-browser.`
  3. Take note of the Jupyter port number that has been assigned (the four digits in the X's in http://127.0.0.1:XXXX) and the token (the long string after `token=`).
  4. Open a new terminal window or tab on your local computer. In this terminal we set up an SSH tunnel.
  5. Pick a random number YYYY between 8000 and 9000. This will be our SSH port number for the tunnel. Try another number if something fails.
  6. On your local machine, type `ssh -A -L YYYY:localhost:XXXX 1234567@gemini.science.uu.nl`
  7. Open a browser on your local computer and go to `localhost:YYYY`, where `YYYY` is your chosen portnumber. When asked for a password/token, use the one that you noted in step 2.

More info can be found here: https://github.com/OceanParcels/UtrechtTeam/wiki/How-to-run-parcels-on-lorenz,-gemini-and-cartesius#gemini

</details>




1. You can download the files in this repository by running `git clone https://github.com/daanreijnders/trapping_internal_tides.git`
1. Navigate to the `model` directory: `cd trapping_internal_tides/model/`
1. Build the model with the `make` command



## Model components

You will find the following files in the model directory:

- `*.f` files. These are the source codes for the model. The main code is `micom-ssd.f`, and the other files serve as subroutines.
- `makefile` that contains instructions for the compiler
- `*.o` files that are created after compiling the code with the makefile
- `micom1.x` file created after compiling the code with the makefile. This is executable for running the model
- `*.h` (3x) files where some variables are defined common to various subroutines and the main part
- `micom.in` file setting a few free parameters.
- `thetas` file describing the densities of the 43 layers.



## Running the model

With the command `./micom1.x`the model starts running. **If you are working on the cluster, do not execute this command directly from the shell. Instead, use the queue system. See the *⚠️ Instructions for using the Gemini cluster* above.

Three types of output files will be created;

- 16 `analyse_****` files

- 3 `average_****` files

- 1 `restart_****` file

 The MICOM model is a so-called isopycnic model. It consists of layers with constant density but the depth of the interface between the layers is variable (in contrast to so-called constant z-models, where the density of layers is variable, but their interfaces are fixed in time).

- The analyse-files output the time varying depth of the 43 interfaces at 16 moments in time during the last day (1.5 hour output). This is enough to resolve the dominant frequency of internal waves forced by an external wave with fixed frequency of 12 hrs.
- The average-files describe the 1-day average of the model-state. You only need the last day.
- The restart-file is needed when you want to prolong your run



## Reading and analyzing the model output

An analysis.py is available to read the output files in a python-script and to calculate amplitude and phase of the dominant internal wave. Amplitudes must be divided by <img src="https://render.githubusercontent.com/render/math?math={9.806\times10^4}"> to convert them to meters (<img src="https://render.githubusercontent.com/render/math?math={g = 9.806 \ \rm{m/s^2}}"> is the gravitational acceleration); the phase is defined between <img src="https://render.githubusercontent.com/render/math?math={-\pi}"> and <img src="https://render.githubusercontent.com/render/math?math={\pi}">.



## Assignment

In the first practical you will try to recover some of the figures in Drijfhout and Maas (2007) and possibly other figures from the same 4 runs. We ask you to

- Motivate your choice of figures (choose 6-8).

- Explain why you must adapt <img src="https://render.githubusercontent.com/render/math?math={N}"> and total runtime when <img src="https://render.githubusercontent.com/render/math?math={f}"> changes
- Why you see trapping of internal waves in the present set-up and not in other set-ups (different <img src="https://render.githubusercontent.com/render/math?math={f}"> values and bottom profiles)
- What is the story you want to tell (see the first bullet) and what are your main conclusions?

 In practical 1 you do not need to change the set-up (but you may and can win extra brownie points if you do). In practical 2 you will have to change the set-up. Here are the places you will have to make changes in the code (although the code is in fortran the required changes are so simple that you can make them without detailed understanding of the code).

- To change the bottom profile, in `cyclo.f` you must edit `pmer` (line 57) and `poos` (line 49) to define the bottom profile in respectively length and width.
- To change <img src="https://render.githubusercontent.com/render/math?math={f}"> you must edit `geopar.f` (line 22-23)
- To change <img src="https://render.githubusercontent.com/render/math?math={N}"> you must edit `thetas`. The values in thetas represent the potential densities, and are determined as (<img src="https://render.githubusercontent.com/render/math?math={\sigma_0)-1000/1000}">. That means that <img src="https://render.githubusercontent.com/render/math?math={1000*\theta}"> runs from 26.0 till 28.071 <img src="https://render.githubusercontent.com/render/math?math={\rm{kg}/{m}^3}">. 
- To change the length of the run you must edit `micom.in`. The 5 values are explained in `micom_ssd.f` where `micom.in` is read on line 27. You can search for their names to understand what they steer. The first 2 values refer to `day1` and `day2` and the model runs from `day1` to `day2`.

Note that you may need to recompile the program using `make`.

In practical 2 we ask you to configure at least 1 of the 3 other set-ups discussed in Drijfhout and Maas (2007) (if you have already done so in practical 1 we ask you to choose a second one), and at least 1 configuration that was **NOT** discussed in Drijfhout and Maas (2007). You may think of the following options: 

- What happens when you change the forcing frequency? Are the results invariant for this parameter? And what happens if you add the new and old frequency in the forcing? How linear is the response?
  You can change the forcing frequency in `boundpb.f` and `boundvb.f` on line 30 (`sin(tsec*2.*pi/(12.*3600.))`).
- What happens if you choose steeper or less steep bottom profiles in the cross-channel direction? What happens if you add random small-scale perturbations to the bottom?
- What happens if you turn the channel, or to make it easy, make <img src="https://render.githubusercontent.com/render/math?math={f}"> a variable function of *x*?
- What happens if you make <img src="https://render.githubusercontent.com/render/math?math={N}"> stronger or weaker or no longer a constant but a function of *z*?
- Maybe you want to change another parameter. Which one would you choose (try to argue, even if you have no time to do this)?

Again we ask you:

- To motivate your choice of figures (choose 2-4) for each of the 2 runs with 1 parameter changed.

- What is the story you want to tell (see bullet 1) and what are your main conclusions?