#!/bin/sh
rsync --progress -auvz --delete ../list/    mucat:/home/mucat/e_202302data/SDDdata/list/
rsync --progress -auvz --delete ../root/    mucat:/home/mucat/e_202302data/SDDdata/root/
