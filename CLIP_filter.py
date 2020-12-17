import sys
import argparse



def main():
    parser = argparse.ArgumentParser(description='finds gga within CLIP-SEQ')
    parser.add_argument('clip', type=str, help='Control <.react> file')
    parser.add_argument('gga',type=str, help='Experimental <.react> file')

    args = parser.parse_args()

    CLIP_GGA = []

    with open(args.clip,"rU") as clip:
        firstline_clip = clip.readline()
        with open(args.gga,'rU') as gga:
            firstline_gga = gga.readline()

            for line in gga:
                line = line.strip()
                line = line.split(",")

                transcript_gga = str(line[0])
                if line[1] == 'NA':
                    FPUTR = line[1]
                else:
                    FPUTR = int(line[1])
                CDS = int(line[2])
                start = int(line[4])
                end = int(line[5])
                length_GGA = len(line[6])

                if FPUTR == 'NA':
                    gga_dist_from_TSS = start
                elif FPUTR > start:
                    gga_dist_from_TSS = start - FPUTR
                else:
                    gga_dist_from_TSS = start - FPUTR

                for item in clip:
                    item = item.strip()
                    item = item.split(",")

                    transcript_CLIP = str(item[1])
                    clip_dist_from_TSS = item[5]
                    length_CLIP = len(item[4])

                    if item[5] != 'NA':
                        clip_dist_from_TSS = int(item[5])

                        if clip_dist_from_TSS > 0 and gga_dist_from_TSS > 0:
                            distance = clip_dist_from_TSS-gga_dist_from_TSS

                            if transcript_gga == transcript_CLIP:
                                if length_CLIP > distance+12:
                                    if ",".join(line) not in CLIP_GGA:
                                        CLIP_GGA.append(",".join(line))


                        elif clip_dist_from_TSS < 0 and gga_dist_from_TSS > 0:
                            distance = abs(clip_dist_from_TSS)+gga_dist_from_TSS

                            if transcript_gga == transcript_CLIP:
                                if length_CLIP > distance+12:
                                    if ",".join(line) not in CLIP_GGA:
                                        CLIP_GGA.append(",".join(line))


                        elif clip_dist_from_TSS < 0 and gga_dist_from_TSS < 0:
                            distance = abs(clip_dist_from_TSS)-abs(gga_dist_from_TSS)

                            if transcript_gga == transcript_CLIP:
                                if length_CLIP > distance+12:
                                    if ",".join(line) not in CLIP_GGA:
                                        CLIP_GGA.append(",".join(line))


                        elif clip_dist_from_TSS == 0 and gga_dist_from_TSS == 0:

                            if transcript_gga == transcript_CLIP:
                                if length_CLIP > 12:
                                    if ",".join(line) not in CLIP_GGA:
                                        CLIP_GGA.append(",".join(line))

                clip.seek(0)
                firstline_clip = clip.readline()


    NO_GGA = []

    with open(args.gga,'rU') as gga:
        firstline_gga = gga.readline()
        for line in gga:
            line = line.strip()

            if line not in CLIP_GGA:
                NO_GGA.append(line)






    with open("CLIP_GGA_foot.csv","w") as outp:
        outp.write(firstline_gga)
        outp.write("\n")
        for item in CLIP_GGA:
            outp.write(item)
            outp.write("\n")

    with open("NON_GGA_foot.csv","w") as outp:
        outp.write(firstline_gga)
        outp.write("\n")
        for item in NO_GGA:
            outp.write(item)
            outp.write("\n")





if __name__ == '__main__':
    main()
