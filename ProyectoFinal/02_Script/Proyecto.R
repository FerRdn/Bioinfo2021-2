####################
##### Proyecto #####
####################

### Análisis bioinformático con R base.

## Comparar secuencias pero usa Biostrings para leer las secuencias

seq1 <- "ATGA"
seq2 <- "TTCT"
# CGTACGTA
# ATGA tiene codón de paro e inicio jaja

class (seq1) # Caracteres

## 1. Empezar con una secuencia corta de nucleótidos.
# Pide dos secuencias cortas a ingresar, ej "ATCG"

seq1 <- readline (prompt = "Ingresa una secuencia corta: ")
seq2 <- readline (prompt = "Ingresa otra secuencia corta: ")


## 2. Largo de la secuencia (número de nc).

print (paste ("La secuencia 1 tiene", nchar (seq1), "nucleótidos."))
print (paste ("La secuencia 2 tiene", nchar (seq2), "nucleótidos."))

# Antes usaba:
# > seq1c <- as.character (seq1)
# > seq1c <- strsplit (seq1c, "") [[1]]
# Para separar los caracteres y contarlos como nc, pero
# > nchar cuenta los caracteres del mismo vector sin separarlos.


## 3. Comparación de tamaños y diferencias de tamaños
# Si ambas secuencias son del mismo tamaño, indica cuántos nc tienen.
# Si la secuencia 1 es menor, se resta la más grande (2) de la más corta (1)
# Para que de la diferencia de nc entre las secuencias.
# Si la secuencia 2 es menor, se resta la más grande (1) de la más corta (2)
# Para que de la diferencia de nc entre las secuencias.

if (nchar (seq1) == nchar (seq2)) {
  print (paste ("Ambas secuencias tienen", nchar (seq1), "nucleótidos."))
} else if (nchar (seq1) < nchar (seq2)) {
  print (paste ("Las secuencias difieren por", nchar (seq2) - nchar (seq1), "nucleótido(s)."))
} else {
  print (paste ("Las secuencias difieren por", nchar (seq1) - nchar (seq2), "nucleótido(s)."))
}


## 4. Mismatches, similaridad o igualdad
## En secuencias del MISMO tamaño
# Si las dos secuencias son exactamente iguales, pues son iguales y ya.
# Si difieren, se indica cuáles bases son las que difieren (no cuenta total)
# Y muestra dónde se encuentran con seq1 == seq2.

if (seq1 == seq2) {
  print ("Las dos secuencias son iguales.")
} else {
  seq1c <- as.character (seq1)
  seq1c <- strsplit (seq1c, "") [[1]]
  seq2c <- as.character (seq2)
  seq2c <- strsplit (seq2c, "") [[1]]
  compara <- seq1c == seq2c
  print (paste ("Las secuencias difieren en", adist (seq1, seq2), "nucleótido(s)."))
  print (compara)
}

# Aquí, si quieres, indica exactamente en cuales nc no coinciden.
respuesta <- 1
respuesta <- readline (prompt = "¿Quieres saber en cuáles difieren?
                         Sí = 1; No = 0")
if (respuesta == 1) {
  print (paste ("Las secuencias difieren en el nucleótido", grep ("FALSE", compara)))
} else {
  print ("Pues te lo pierdes.")
}


## 5. Contar las cuatro bases y obtener su porcentaje
# Ciclo for porque el límite son los caracteres en los vectores
# Contar cuántas A, T, C o Gs hay en los vectores
# E imprimir ese numerito
# Luego obtener el porcentaje con ese numerito / length de seq x 100.

adenina <- 0
timina <- 0
citocina <- 0
guanina <- 0


for (adenina in seq1) {
  if (seq1[1] == nchar("A")) {
  adenina <- adenina + 1
  print (adenina)
  }
}


seq1c [2] # Objeto dos del vector = "T"

length (setdiff (seq1, seq2))

countPattern ("NA", match (seq2c, seq1c))

seq1c == seq2c

grep ("FALSE", compara) # Busca, no cuenta

seq1c [1]

adist (seq1, seq2)

char (seq1) %in% char (seq2)
seq2c %in% seq1c # Toma a todo el vector como un solo caracter

str_count ("NA", match (seq1c, seq2c))

compara <- compareStrings (seq1c, seq2c)
countPattern ("?", compara)

## 6. Contar los codones de inicio y de paro.
# No separa por codones, entonces como en el caso de seq1
# se tiene tanto un codón de inicio como de paro
# cuando sólo debería de haber de inicio porque está primero.

# Inicio

inicio <- function (seq1) {
  seq1cod <- strsplit (seq1, "") [[3]]
  grep ("ATG", seq1cod)
}

inicio (seq1cod)
inicio (seq2cod)

# Paro

paro <- function (seq1) {
  seq1cod <- strsplit (seq1, "") [[3]]
  grep ("TAA", seq1cod)
  grep ("TAG", seq1cod)
  grep ("TGA", seq1cod)
}

paro (seq1)
paro (seq2)


## 7. Obtener la cadena de RNA
# La T se sustituye por una U en las secuencias.
# Con la función aplica para cualquier secuencia.

RNA <- function (seq1) {
  gsub ("T", "U", seq1)
}

RNA (seq1)
RNA (seq2)

seq1RNA <- gsub ("T", "U", seq1)
seq1RNA

seq2RNA <- gsub ("T", "U", seq2)
seq2RNA


## 8. Obtener la cadena complementaria
# Para A es T, para T es A, para C es G y para G es C.
# Se aplica la misma de arriba.

Comp <- function (seq1) {
  gsub ("A", "T", seq1)
  gsub ("T", "A", seq1)
  gsub ("C", "G", seq1)
  gsub ("G", "C", seq1)
}
# Está re sustituyendo el cambio y por eso se queda igual
# Tal vez algún ciclo con el cual cambie base por base?
Comp (seq1)
Comp (seq2)


## 9. Obtener la cadena reverso
# La primera base pasa a ser la última,
# La segunda pasa a ser la penúltima, y así.




## 10. Obtener la cadena de aá
# Separar por codones
# asignar valores a los codones pero que siga formando un vector
# Y que sea en orden.


## 11. Sacar la cantidad de aá

nchar ()

## 12. Comparar tamaños de las cadenas de aá

