###################################
##### ANÁLISIS BIOINFORMÁTICO #####
###################################

##### Proteína KRT9

## Librerías
library (Biostrings) # Secuencias
library (ggplot2) # Algunas gráficas
library (msa) # Alineamiento
library (ape) # Árboles
library (seqinr) # Matrices
library (ggtree) # Árboles más bonitos


## Abrir el archivo concatenado
KRT9 <- readAAStringSet (choose.files ()) # Archivo KRT9.faa
KRT9


## Cambiar por nombres comunes
names (KRT9)
names (KRT9) <- c ("Panda Gigante", "Dingo", "Perro", "Nutria Marina",
                   "León Marino", "Gato", "Humano", "Hiena", "Foca",
                   "Nutria Común", "Lince de Canadá", "Lince Ibérico",
                   "Armiño", "Hurón", "León", "Tigre", "Gato leopardo",
                   "Suricata", "Oso gris", "Oso Polar", "Alpaca",
                   "Zorro Ártico", "Zorro común", "Lobo Marino")
KRT9


## Tamaño de las secuencias
width (KRT9) # mushos númeritos
mean (width (KRT9)) # Media de 648.125 aá


## Frecuencia de aá
KRT9aa <- alphabetFrequency (KRT9) # Ver los aá que contiene
KRT9aa # Secuencias ricas en Gly
min (KRT9aa [ , 8]) # Mín 75 aá de Gly
max (KRT9aa [ , 8]) # Máx 269 aá de Gly

KRT9aa1 <- as.data.frame (KRT9aa) # Data frame para hacer gráfica

especies <- c ("Panda Gigante", "Dingo", "Perro", "Nutria Marina",
               "León Marino", "Gato", "Humano", "Hiena", "Foca",
               "Nutria Común", "Lince de Canadá", "Lince Ibérico",
               "Armiño", "Hurón", "León", "Tigre", "Gato leopardo",
               "Suricata", "Oso gris", "Oso Polar", "Alpaca",
               "Zorro Ártico", "Zorro común", "Lobo Marino")
KRT9aa1 <- cbind (KRT9aa1, especies)
# Crear vector con los nombres comunes y agregarlos como columna a los aá

KRT9aa1

pdf ("03_Plots/GlyEnKRT9.pdf")
frecaa <- ggplot (KRT9aa1, aes (x = especies, y = G, fill = especies)) +
  geom_bar (stat = "identity") +
  xlab ("Especies") +
  ylab ("Glicinas") +
  ggtitle ("Contenido de Glicina") +
  theme (axis.text.x = element_text (angle = 90))
frecaa
dev.off ()


## Alineamientos Múltiples
# Con MUSCLE
KRT9AM <- msa (KRT9, "Muscle")
print (KRT9AM, show = "alignment")

# Con ClustalOmega
KRT9AO <- msa (KRT9, "ClustalOmega")
print (KRT9AO, show = "alignment")

# Con ClustalW
KRT9AW <- msa (KRT9, "ClustalW")
print (KRT9AW, show = "alignment")


## Árboles Filogenéticos
# Con MUSCLE
KRT9M <- msaConvert(KRT9AM, type="seqinr::alignment") # Convierte la alineación para hacer el árbol
KRT9M1 <- dist.alignment (KRT9M, "identity") # Alineamiento para matriz de distancias ¿?
as.matrix (KRT9M1) # Matriz bien estructurada
KRT9M2 <- nj (KRT9M1) # Contrucción del árbol con Neighbour Joining
plot (KRT9M2, main = "KRT9 con Muscle") # Visualización

# De ClustalOmega
KRT9O <- msaConvert(KRT9AO, type="seqinr::alignment")
KRT9O1 <- dist.alignment (KRT9O, "identity")
as.matrix (KRT9O1) 
KRT9O2 <- nj (KRT9O1)
plot (KRT9O2, main = "KRT9 con ClustalOmega")

# De ClustarW
KRT9W <- msaConvert(KRT9AW, type="seqinr::alignment")
KRT9W1 <- dist.alignment (KRT9W, "identity") 
as.matrix (KRT9W1) 
KRT9W2 <- nj (KRT9W1) 
plot (KRT9W2, main = "KRT9 con ClustalW")


# Árboles Filogenéticos Avanzados
# Con MUSCLE
pdf ("03_Plots/KRT9_MUSCLE.pdf")
KRT9MT <- ggtree (KRT9M2, aes (branch.length = "Scale", color = branch.length)) + # Asigna colores a ramas
  scale_color_continuous (name ='Scale', limits = c (0, 1.5), oob = scales::squish,
                         low = 'deeppink4', high = 'lightpink') + # Le da color difundido a la escala 
  geom_tiplab (size = 4, color = "deeppink4") + # Texto de las ramas de tamaño y color
  geom_point2 (aes (subset = !isTip), shape = 21, size = 1.5, fill = 'deeppink4') + # En todos los nodos pone figuras
  theme_tree2 ("white") # Fondo blanco
KRT9MT
dev.off ()

# Con ClustalOmega
pdf ("03_Plots/KRT9_ClustalO.pdf")
KRT9OT <- ggtree (KRT9O2, aes (branch.length = "Scale", color = branch.length)) + # Asigna colores a ramas
  scale_color_continuous (name ='Scale', limits = c (0, 1.5), oob = scales::squish,
                          low = 'deeppink4', high = 'lightpink') + # Le da color difundido a la escala 
  geom_tiplab (size = 4, color = "deeppink4") + # Texto de las ramas de tamaño y color
  geom_point2 (aes (subset = !isTip), shape = 21, size = 1.5, fill = 'deeppink4') + # En todos los nodos pone figuras
  theme_tree2 ("white") # Fondo blanco
KRT9OT
dev.off ()

# Con ClustalW
pdf ("03_Plots/KRT9_ClustalW.pdf")
KRT9WT <- ggtree (KRT9W2, aes (branch.length = "Scale", color = branch.length)) + # Asigna colores a ramas
  scale_color_continuous (name ='Scale', limits = c (0, 1.5), oob = scales::squish,
                          low = 'deeppink4', high = 'lightpink') + # Le da color difundido a la escala 
  geom_tiplab (size = 4, color = "deeppink4") + # Texto de las ramas de tamaño y color
  geom_point2 (aes (subset = !isTip), shape = 21, size = 1.5, fill = 'deeppink4') + # En todos los nodos pone figuras
  theme_tree2 ("white") # Fondo blanco
KRT9WT
dev.off ()
