# LBM_Fluid-Flow_Porous-Media

Project Description:

This project is based on my Master of Science Thesis which is titeld "Pore  Scale  Simulation  of  Wettability  and  Interfacial  Tension  on  Two  Phase  Flow through Porous Media by Phase Field Free Energy Lattice Boltzmann Method". In this project two-phase flow, with defined Wettability and Interfacial Tension (IFT), is considred in 2D porous medium; the proposed 2D porous medium is uploaded as "200_200_Periodic" name in this repository. In addition, for simulating various Wettability-amount an Excel file is uploaded as .... which should be used which the procedure is described in following section. 

How to Run the Model:

For running the project, you should download the repository on your PC.  simply open it via Colab, then running the first cell; it will automatically download the dataset, auto-label them, and finally after pre- processing of images they will be fed into the model and training is performed. (Note: Some of the details related to this program are listed in "Extra Explanation").

You can check the model by running "Model Prediction" which is prepared at the end of program. This code will allow you to choose one or more files from your file system, upload them, and run them through the model, giving an indication of whether the object is a horse or a human. please note that as this project is simply designed with low number of training images, it is not always correct.

Extra Explanation:

    For training the model "Horses or Humans" dataset which is prepared by Laurence Moroney (Instructor of Coursera) is used.
    This dataset contains over a thousand images of horses and humans . Some of the images are showed to get a better sense of what they look like.
    ImageDataGenerator class is used to prepare the dataset, so it can be fed to a convolutional neural network.
    The model is developed through TensorFlow, and totally 5 Convolution layers are used (The model architecture can be checked via model.summary()).
    The model is trained with "binary_crossentropy" loss and "RMSprop" optimizer.
    Preprocessing of images is done by ImageDataGenerator; the images are auto-labeled, and they are fed to the model in baches of 128 and size of 300*300
    Training is performed for 15 epochs
