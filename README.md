# ImageToADINA

#### Author Mario Malic, August 2020

Contact mail: mario@mariomalic.com

Alternative contact mail: mario.malic@gmail.com

#### Credits:
- Function 'rgb2hex' written by Chad A. Greene, sourced from MATLAB File Exchange
##### Recommendations:
- There are few options to reduce the time to build the model in addition to environment control settings:  
- Edit - Memory Usage and give ADINA-AUI extra RAM to work with.  
- File - Stop Recording if it is turned on
- Hide message window on the bottom of AUI window
- On my machine 161 x 240 image takes 7 minutes, 402 x 600 takes an hour 515 x 768 six hours.
- ADINA will most likely freeze (not respond) during the model buildup, have patience.
- For the larger band tables 5MB+, it is faster to copy and paste rather to import it

#### Running the model
- Download rgb2hex.m from MATLAB File Exchange
- Put ImageToADINA.m, rgb2hex.m and desired image in the same folder
- Run the code
- Choose the image you would like to create in ADINA when prompted
- Output is the ADINA model input file (.in) and custom band table (.txt)
- Import the input file and solve it
#### Postprocessing the model
- Load the porthole file and create a band plot of axial strain values
- Import Band_File_Name.txt as custom band table for strain values
