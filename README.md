# LBM_Fluid-Flow_Porous-Media

Project Description:

This project is based on my Master of Science Thesis which is titeld "Pore  Scale  Simulation  of  Wettability  and  Interfacial  Tension  on  Two  Phase  Flow through Porous Media by Phase Field Free Energy Lattice Boltzmann Method". In this project two-phase flow, with defined Wettability and Interfacial Tension (IFT), is considred in 2D porous medium; the proposed 2D porous medium is uploaded as "200_200_Periodic" name in this repository. In addition, for simulating various Wettability-amount an Excel file is uploaded as .... which should be used which the procedure is described in following section. 

How to Run the Model:
For running the project, you should follow the following instructions:
1. Download the repository on your PC. 
2. Then please copy and paste the directory of defined porous medium "200_200_Periodic.bmp" on your PC in line number 38 of "Relative_Permeability_Porous_Media.m" code (It's already "C:\Users\SONY\Desktop\LBM Study\Porous Media\200_200_Periodic.bmp"). 
3. For changing the Wettability value to desired one you can change the coefficient of formula in line number 98 of code based on column "D" of uploaded excel file, named "Wetability" for any desired value of contact angle (The coefficient is already set as 0.5863 for contact angle of 30 degrees). 
4. For setting the saturation value of Wet Phase you can change it in line number 17 of the code, and the Non-Wet saturation will be updated automatically. 
5. Then you can run the code.


Extra Explanation:

    For training the model "Horses or Humans" dataset which is prepared by Laurence Moroney (Instructor of Coursera) is used.
    This dataset contains over a thousand images of horses and humans . Some of the images are showed to get a better sense of what they look like.
    ImageDataGenerator class is used to prepare the dataset, so it can be fed to a convolutional neural network.
    The model is developed through TensorFlow, and totally 5 Convolution layers are used (The model architecture can be checked via model.summary()).
    The model is trained with "binary_crossentropy" loss and "RMSprop" optimizer.
    Preprocessing of images is done by ImageDataGenerator; the images are auto-labeled, and they are fed to the model in baches of 128 and size of 300*300
    Training is performed for 15 epochs
