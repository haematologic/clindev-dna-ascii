from collections import defaultdict, Counter

def loadfile(path):
    '''open file for reading as buffered text stream and return iterable lines object
    allows potentially unlimited input file size by only reading line-wise'''
    try:
        with open(path, 'rt') as f:
            lines = list(f)
    except IOError:
        raise
    return lines

def strip_line(line):
    '''strip line breaks and return string of chars as a list'''
    return list([''.join(x.strip('\n') for x in line)][0])

def bin_map(seq):
    '''map DNA to binary per given keys'''
    mappings = {'A':'0','C':'0','G':'1','T':'1'}
    return list(map(lambda x: mappings.get(x) , seq))

def parse_dna_lines(line):
    '''manually handle FASTA format returning only lines not starting ">" and map DNA to binary'''
    stripped_line = strip_line(line)
    if stripped_line[0] == '>':
        return
    else:
        return bin_map(stripped_line)

def parse_lines(lines):
    '''generator yielding binary strings read from seq file stream iterator'''
    for line in lines:
        bin_line = parse_dna_lines(line)
        if bin_line:
            yield (bin_line)

            
def buffer_lines(buffer, lines):
    '''fill buffer with line of binary until end of file (EOF)'''
    for line in lines:
        binline = line
        buffer += binline
        break
    return buffer
        

def parse_bytes_by_offset(lines):
    '''main function for parsing sequental groups of 8 bits, offset by 0-7 bits inclusive'''
    # initialise defaultdict – offers some useful enhancements over standard dictionary type
    # e.g. creation of keys if not existing when called via append method as below
    all_offsets = defaultdict(list)  

    lines = parse_lines(lines)
    
    # initialise empty buffer with line of binary string parsed from file'''
    buffer = buffer_lines([], lines)
    bits_parsed = 0

    # now loop until buffer empty
    while buffer:
        # get line of binary from generator if less than 16 bits to parse'''
        if len(buffer)<16:
            buffer = buffer_lines(buffer, lines)
            if not buffer:
                break
        
        # parse bytes with start bit offset by by 0 to 7
        for offset in range(0,8):
            # calculate global_offset from absolute start of binary string and current offset / bits_parsed
            global_offset = bits_parsed + offset
            
            # slice buffer to obtain 8 bits starting at specified offset
            bin_byte_str = ''.join(buffer[offset:offset+8])
            
            # calculate decimal integer equivalent and actual charater
            dec_char = int(bin_byte_str, 2)
            char = chr(dec_char)

            # restrict 'hits' to limited range of ASCII chars
            if 32 <= dec_char <= 127 or dec_char in {11, 12, 15}:
                # append tuple containing global_offset and character
                offset_char = (global_offset, char)
                all_offsets[offset].append(offset_char)  
        
        # slice buffer to move along by 8 bits
        buffer = buffer[8:]
        
        # update position of current bits_parsed
        bits_parsed += 8
    
    return all_offsets  # return collections.defaultdict


def calc_char_freqs(all_offsets):
    '''do frequency analysis to find offset containing most "hits:'''
    ascii_char_freqs = []
    
    # loops through offsets 0 to 7 in parsed_bytes
    for offset in all_offsets:
        
        # retrive characters from tuple containing global offset and character
        chars = [tup[1] for tup in all_offsets[offset]]
        
        # initialise collections.Counter datatype to aid freq analysis
        c = Counter(chars)
        
        result = offset, len(all_offsets[offset]), c.most_common(n=5)
        ascii_char_freqs.append(result)
    
    return ascii_char_freqs


def main(seq_file = './data/id100228968.seq'):
    lines = loadfile(seq_file)
    
    # call main function to parse bytes line-wise
    all_offsets = parse_bytes_by_offset(lines)
    
    # calc frequencies of ascii chars
    ascii_char_freq = calc_char_freqs(all_offsets)

    # sort offsets by most hits
    ascii_char_freq_sorted = sorted(ascii_char_freq, key=lambda x: x[1], reverse=True)
    print('Frequency analysis of bit offsets')
    [print(x, end='\n') for x in ascii_char_freq_sorted]
    print()
    
    print('Results of offset having most hits')
    for tup in all_offsets[ascii_char_freq_sorted[0][0]]:
        print(tup[1], end='')

main()