import xlrd
import random
import collections

rows,seq,pos = [],[],[]

with open('intrinsic.csv','r') as inp:
    for line in inp:
        pos.append(str(line.split(",")[0]))

# accessing Column 'C' in this example
COL_IDX = 1

book = xlrd.open_workbook('intrinsic.xls', formatting_info=True)
first_sheet = book.sheet_by_index(0)

for row_idx in range(first_sheet.nrows):
  text_cell = first_sheet.cell_value(row_idx, COL_IDX)
  text_cell_xf = book.xf_list[first_sheet.cell_xf_index(row_idx, COL_IDX)]

  # skip rows where cell is empty
  if not text_cell:
    continue
  #print text_cell,

  text_cell_runlist = first_sheet.rich_text_runlist_map.get((row_idx, COL_IDX))
  if text_cell_runlist:
    #print '(cell multi style) SEGMENTS:'
    segments = []
    for segment_idx in range(len(text_cell_runlist)):
      start = text_cell_runlist[segment_idx][0]
      # the last segment starts at given 'start' and ends at the end of the string
      end = None
      if segment_idx != len(text_cell_runlist) - 1:
        end = text_cell_runlist[segment_idx + 1][0]
      segment_text = text_cell[start:end]
      segments.append({
        'text': segment_text,
        'font': book.font_list[text_cell_runlist[segment_idx][1]]
      })
    # segments did not start at beginning, assume cell starts with text styled as the cell
    if text_cell_runlist[0][0] != 0:
      segments.insert(0, {
        'text': text_cell[:text_cell_runlist[0][0]],
        'font': book.font_list[text_cell_xf.font_index]
      })

    i = 0
    for segment in segments:
        #print segment['text'],
        #print "\n"

        if i % 2 != 0:
            rows.append(row_idx)
            seq.append(str(segment['text']))
        i += 1

  else:
    #print '(cell single style)',
    #print 'italic:', book.font_list[text_cell_xf.font_index].italic,
    #print 'bold:', book.font_list[text_cell_xf.font_index].bold
    pass


with open('duplex_len_new.csv','w') as outp:
    q=0
    for item in seq:
        outp.write(str(pos[q])+","+item+"\n")
        q+=1
