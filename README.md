# AVDASI-4
AVDASI 4 (Aerospace Vehicle Design and Systems Integration) – Group Design Project is a 40 Credit Points Design Project set by Airbus Defence & Space. We worked as a group of seven students and we were tasked with designing a satellite to study the composition of Venus, based on the Envision satellite developed by ESA and set to launch in 2032.
The deliverables included two presentations from the group given to academics and Airbus employees, two technical group reports (including a final report of 56 pages), and an individual presentation focused on our subsystem. I received an overall mark of 77.8%, including 89% for the individual presentation.
This project aimed to immerse us in a real industry project to fully understand the daily challenges faced by engineers.


My specific role was the power subsystem. The initial solar panels were sized in MATLAB through the power requirements given by Airbus Defence & Space. Then STK ANSYS was used to model the power produced by the satellite at all times.
To model the power production on ANSYS STK the following steps have been followed:
  1. CAD modelling of the satellite with the solar panels on AutoDesk Inventor
  2. Export as an STD file
  3. Open in Blender to adjust the orientation of the satellite
  4. Export as a .dae file
  5. Open in Visual Studio and modify the code to ensure ANSYS STK can recognise the solar cells and code the articulations of the solar panels
  6. Upload the .dae file to ANSYS STK


This repository includes:

  - Power calculations done on MATLAB. It does not include all the MATLAB codes used, only the two most relevant ones.
  - .dae file uploaded on ANSYS STK to model the satellite (after the modifications explained in step 5 above).
