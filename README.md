# Trapping of Internal Tides in a Channel Model

Two tutorials of the Utrecht University course NS-MO447M Waves in Geophysical Fluids will be devoted to recovering and extending the result published in the paper

> Drijfhout, S., & Maas, L. R. M. (2007). Impact of Channel Geometry and Rotation on the Trapping of Internal Tides. *Journal of Physical Oceanography*, *37*(11), 2740–2763. https://doi.org/10.1175/2007JPO3586.1

which you can download here: https://webspace.science.uu.nl/~maas0131/files/drijfhoutmaas07jpo[smallpdf.com].pdf

You will do so using the simple ocean circulation model named “Micom”, written in Fortran, which is the main computer language for nearly all ocean, atmosphere, and climate models.

This repository holds the model code and files that you can use for analysis. This README also contains the instructions for a small report that forms the fifth assignment of the course

##  Downloading and building the model

You are used to working with Python, which is an interpreted programming language. This means that your instructions are not immediately understood by your computer, but are interpreted and translated to machine code on-the-fly as you execute a script. In contrast, Fortran is a compiled programming language. Before you run a program, your entire code needs to be translated to machine code first. This extra step allows your code to run incredibly efficiently whenever it is executed.

To run the model on your machine, you will first need to compile it. This is called *building* the model. The model should be able to build on MacOS and Linux with the Fortran compiler *gfortran* (part of the *GCC* compiler collection) installed. This is default on Gemini, the Science department’s computer cluster. If you are using MacOS, make sure to have the Command Line Tools installed: open the Terminal and type `xcode-select --install`.

<details>
  <summary>Instructions for using the Gemini cluser</summary>

  ### Logging in
  1. Open a Terminal
  2. Connect to the Gemini cluster by typing `ssh XXXXXXX@gemini.science.uu.nl` using your Solis-ID in place of the XXs
  3. Type your Solis-ID password
  4. You're in!

  ### Downloading files
  1. You can download the files in this repository by running `git clone
  
  ### Submitting a job
  
</details>


On you workstation you have access to a series of files:

·   *.f. files (these are the source codes for the model, a main part, micom-ssd.f, and subroutines)

·   1 makefile

·   *.o files that are created after compiling the code with the makefile

·   1 micom1.x file created after compiling the code with the makefile

·   3 *.h files where some variables are defined common to various subroutines ad main part

·   1 micom.in file setting a few free parameters.

·   1 thetas file describing the densities of the 43 layers.

 

The model is configured as a channel with an open boundary at one, oceanward side, and a continental slope at the opposite side. The length of the channel is 1200km, the width 191.25 km. There are 43 layers of 100-m depth, so that total depth in the middle of the channel is 4300 m. 

The flow (internal tide) is forced at the open boundary by a barotropic, semidiurnal tide with a typical period of 12 h. 

The paper discusses 4 experiments. In two experiments the earth’s rotation is neglected and the Coriolis frequency is set to 0. In the other two experiments we use an *f* plane with *f* = 10-4 s-1. Both cases are run with two channel configurations: in one case the channel has a flat bottom and vertical walls in the cross-channel direction; in the other case the channel has a parabolic cross-channel bottom profile with maximum depth of 4300 m and sloping sidewalls. When f = 0 the buoyancy frequency is chosen as *N* = 3.05 10-3 s-1, when f = 10-4 s-1, *N* = 2.2 10-3 s-1. 

On your workstation the model is configured for the case *f* = 10-4 s-1 and *N* = 2.2 10-3 s-1. 

 



 



 

With the command

./micom1.x

The model starts running

Three types of output files will be created;

16 analyse_**** files

3 average_**** files

1 restart_**** file

 

The Micom model is a so-called isopycnic model. It consists of layers with constant density but the depth of the interface between the layers is variable (in contrast to so-called constant z-models, where the density of layers is variable, but their interfaces are fixed in time).

 

·   The analyse-files output the time varying depth of the 43 interfaces at 16 moments in time during the last day (1.5 hour output). This is enough to resolve the dominant frequency of internal waves forced by an external wave with fixed frequency of 12 hrs.

·   The average-files describe the 1-day average of the model-state. You only need the last day.

·   The restart-file is needed when you want to prolong your run

 

An analysis.py is available to read the output files in a python-script and to calculate amplitude and phase of the dominant internal wave. Amplitudes must be divided by 9.806 104 to convert them to meters (g = 9.806 is the gravity constant in m s-2); the phase is defined between -p and p

 

In the first practical you will try to recover some of the figures in Drijfhout and Maas (2007) and possibly other figures from the same 4 runs. We ask you to

·   Motivate your choice of figures (choose 6-8).

·   Explain why you must adapt N and total runtime when f changes

·   Why you see trapping of internal waves in the present set-up and not in other set-ups (different f values and bottom profiles)

·   What is the story you want to tell (see bullet 1) and what are your main conclusions?

 

 

In practical 1 you do not need to change the set-up (but you may and can win extra brownie points if you do). In practical 2 you will have to change the set-up. Here are the places you will have to make changes in the code (although the code is in fortran the required changes are so simple that you can make them without detailed understanding of the code).

 

·   To change the bottom profile, you must edit cyclo.f. pmer and poos define the bottom profile in respectively length and width.

·   To change f you must edit geopar.f

·   To change N you must edit “thetas”

·   To change the length of the run you must edit micom.in. The 5 values are explained in micom_ssd.f where micom.in is read. You can search for their names in micom.in to understand what they steer. The first 2 values refer to day1 and day2 and the model runs from day1 to day2.

 

In practical 2 we ask you to configure at least 1 of the 3 other set-ups discussed in Drijfhout and Mass (2007) (if you have already done so in practical 1 we ask you to choose a second one), and at least 1 configuration that was NOT discussed in Drijfhout and Maas (2007). You may think of the following options:

 

·   What happens when you change the forcing frequency? Are the results invariant for this parameter? And what happens if you add the new and old frequency in the forcing? How linear is the response?

·   What happens if you choose steeper or less steep bottom profiles in the cross-channel direction? What happens if you add random small-scale perturbations to the bottom?

·   What happens if you turn the channel, or to make it easy, make f a variable function of x?

·   What happens if you make N stronger or weaker or no longer a constant but a function of z?

·   Maybe you want to change another parameter. Which one would you choose (try to argue, even if you have no time to do this)?

 

Again we ask you 

·   Motivate your choice of figures (choose 2-4) for each of the 2 runs with 1 parameter changed.

·   What is the story you want to tell (see bullet 1) and what are your main conclusions?