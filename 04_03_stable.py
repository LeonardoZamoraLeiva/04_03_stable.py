#!/usr/bin/env python
# coding: utf-8

import os
import string
import io

ALPHA = string.ascii_letters

# Leer GBks y TXTs en carpetas y subcarpetas
strains_bgcs = []
strains_gbks = []

for entries in os.listdir('./test'):
    if entries.endswith('.txt'):
        strains_bgcs.append(entries)
    elif entries.endswith('.gbk'):
        strains_gbks.append(entries)


def replace_on_line():
    # linea a modificar
    line_replace = content[count].replace('\t', ',')
    # nombre del gen en la linea

    strain_name = strains_bgcs[j].split('.', 1)[0]
    strain_name = strain_name.split('_', 1)[0]

    bgc_name = line_replace.split(',')[0]
    start = line_replace.split(',')[1]
    end = line_replace.split(',')[2]
    direction = line_replace.split(',')[3]

    # reeplaazar el nombre del gen en la linea por el nombre del BGC
    line_replace = strain_name + ',' + bgc_name + ',' + start + ',' + end + ',' + direction

    # Reemplazar '+' y '-' por '1' y -'1', respectivamente
    if '+' in line_replace.split(',')[4]:
        line_replace = line_replace.replace('+', '1,gene{},{},{}\n', 1).format(gen_en_BGC, j + 1,
                                                                               smcogs(lista_de_genes))
        lista_de_genes.pop(0)
    elif '-' in line_replace.split(',')[4]:
        line_replace = line_replace.replace('-', '-1,gene{},{},{}\n', 1).format(gen_en_BGC, j + 1,
                                                                                smcogs(lista_de_genes))
        lista_de_genes.pop(0)

        # Escribir el archivo
    out_file.write(line_replace)
    print(line_replace)


def smcogs(list_of_genes):
    # clasificaciones antiSMASH
    sm_cogs = ['biosynthetic', 'regulatory', 'transport', 'other']

    # variable out debe ser global para poder escribirlo fuera

    global gen_anterior

    # strip lines del archivo .gbk con la anotacion
    for m in range(len(strains_gbks)):
        if strains_bgcs[j].split('.', 1)[0] in strains_gbks[m].split('.', 1) or strains_bgcs[j].split('_', 1)[0] in \
                strains_gbks[m].split('.', 1):
            with io.open('./test/{}'.format(strains_gbks[m]), 'r', encoding='utf-8') as strain:
                gbk_lines = strain.readlines()

            gbk_lines = [b.strip() for b in gbk_lines]

            # para cada linea del archivo

            gene_first = False
            product_first = False

            for p in range(len(gbk_lines)):
                linea_revisada = p
                reference_count = 0
                while not product_first and not gene_first:
                    if 'REFERENCE' in gbk_lines[linea_revisada]:
                        reference_count =+ 1

                    if '/locus_tag' in gbk_lines[linea_revisada]:
                        if reference_count == 0:
                            pass
                        else:
                            gene_first = True

                    elif 'product' in gbk_lines[linea_revisada]:
                        product_first = True

                    linea_revisada += 1

                break

            for line in range(len(gbk_lines)):

                # salir del ciclo si no quedan mas genes por buscar
                if len(list_of_genes) < 1:
                    break

                # si encontramos el nombre del BGC en la primera linea

                gen_actual_en_comillas = '"' + list_of_genes[0] + '"'
                if len(list_of_genes) > 1:
                    gen_siguiente_en_comillas = '"' + list_of_genes[1] + '"'
                else:
                    gen_siguiente_en_comillas = False

                if gen_actual_en_comillas in gbk_lines[line]:

                    # para cada miembro de la clasificacion sm_cogs
                    tipo_encontrado = False
                    for cog_member in range(len(sm_cogs)):
                        # Esta variable permite no buscar dentro del gen siguiente

                        cds_count = 0
                        tag_count = 0

                        if product_first:
                             contador_para_linea = -1
                        else:
                             contador_para_linea = 1

                        function_to_search = '/gene_functions="' + sm_cogs[cog_member]
                        kind_to_search = '/gene_kind="' + sm_cogs[cog_member]

                        linea_en_revision = line + contador_para_linea

                        while cds_count < 2 and tag_count < 2:
                            # Para salir del ciclo si llegamos al gen siguiente

                            if 'CDS ' in gbk_lines[linea_en_revision]:
                                cds_count += 1

                            if gbk_lines[linea_en_revision] == gbk_lines[line] or not gen_siguiente_en_comillas:
                                pass
                            elif gbk_lines[line].split('=')[0] in gbk_lines[linea_en_revision]:
                                tag_count += 1
                            elif gen_siguiente_en_comillas in gbk_lines[linea_en_revision]:
                                tag_count += 1

                            if 'FEATURES' in gbk_lines[linea_en_revision] or 'protein_id' \
                                    in gbk_lines[linea_en_revision]:
                                if product_first:
                                    linea_en_revision = line + 2
                                    contador_para_linea = 1
                                    tag_count += 1
                                else:
                                    linea_en_revision = line
                                    contador_para_linea = -1
                                    tag_count += 1

                            if 'ORIGIN' in gbk_lines[linea_en_revision]:
                                if product_first:
                                    break
                                else:
                                    contador_para_linea = -1

                            if tag_count > 0:
                                if product_first:
                                    contador_para_linea = 1
                                else:
                                    contador_para_linea = -1

                            # Si encuentra el smcog en la linea, imprimir el smcog correspondiente y guardarlo en out.
                            # Salir del ciclo para seguir con el siguiente gen

                            if kind_to_search in gbk_lines[linea_en_revision] or \
                                    function_to_search in gbk_lines[linea_en_revision]:
                                tipo_encontrado = True
                                break

                            linea_en_revision = linea_en_revision + contador_para_linea

                        # Si encontramos el smcog, salir del ciclo for para seguir con el siguiente gen
                        if tipo_encontrado:
                            gen_anterior = gen_actual_en_comillas
                            if 'biosynthetic' in gbk_lines[linea_en_revision]:
                                if 'biosynthetic-additional' in gbk_lines[linea_en_revision]:
                                    out = 'biosynthetic-additional'
                                else:
                                    out = 'biosynthetic'
                            else:
                                out = sm_cogs[cog_member]

                            return out

                    # Si no encontramos el smcog en ninguna parte, salir del ciclo, asignarlo como other y continuar
                    # con el

                    if not tipo_encontrado:
                        # print('opcion 2')
                        out = 'other'
                        gen_anterior = gen_actual_en_comillas
                        return out


# Crear un archivo en el que se iran escribiendo las lineas, y agregar como header los titulos por linea
working_file = './out/actinosynnema_and_350.txt'
out_file = open('{}'.format(working_file), "a")
out_file.write('strain,antismash_gene,start,end,direction,gene_number,molecule_number,antismash_smcog' + '\n')

# leer todas las lineas y devolerlas sin los caracteres de inicio y final

for j in range(len(strains_bgcs)):
    with io.open('./test/{}'.format(strains_bgcs[j]), 'r', encoding='utf-8') as f:
        content = f.readlines()

    content = [x.strip() for x in content]

    # Buscar los datos que interesan
    for i in range(len(content)):
        if 'Table of genes' in content[i]:
            count = i + 1
            conteo_de_genes = i + 1
            gen_en_BGC = 1
            lista_de_genes = []
            lista_de_genes_2 = []
            gen_anterior = False

            # hacer una lista de los genes para utilizarlos en la funcion
            while content[conteo_de_genes].startswith(tuple(ALPHA)):
                gene_name = content[conteo_de_genes].split('\t', 1)[0]
                lista_de_genes.append(gene_name)

                conteo_de_genes += 1

                if not content[conteo_de_genes] != '':
                    break

            # Copiar toda la tabla con las caracteristicas de los genes
            lista_de_genes_2.append(lista_de_genes[0])
            primer_gen = '"' + lista_de_genes_2[0] + '"'

            while content[count].startswith(tuple(ALPHA)):
                # Funcion para incorporar la linea
                # out = ''

                replace_on_line()
                count += 1
                gen_en_BGC += 1

                # Si es una linea blanca, volver al inicio y reiniciar el conteo de genes a 1 para el siguiente cluster
                if not content[count + 1] != '':
                    gene_en_BGC = 1
                    # break
            break


out_file.close()

# renombrar el archivo a .csv para abrirlo en R
os.rename(r'{}'.format(working_file), r'{}'.format(working_file[0:-4]+'.csv'))

# Ejecutar el script de R para dibujar los genes
os.system('Rscript R_script_prueba.R')
