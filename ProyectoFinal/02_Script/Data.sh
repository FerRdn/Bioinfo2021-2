#!/bin/bash

###### Proyecto Final
##### Fernanda Rendón

# Aquí sólo concatenaré las secuencias de aminoácidos

# En la carpeta de todo el proyecto hay otra carpeta 01_Data
cd ..
ls

# Ahí están las secuencias
cd 01_Data

# Concatenar
cat *faa > KRT9.faa

# Para comprobar
ls

# Listo :)
