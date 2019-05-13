import sys
import re

def parse(fp):
    for line in fp:
        line = line.strip()
        if line[0] == '@':
            continue
            
        dico = parse_columns(line)

        yield dico


def parse_cigar(cigar):
    # Split the cigar line
    errors = re.split(r'([0-9]+[M,D,I])', cigar)
    errors = [val for val in errors if len(val) > 0]

    # Count each error
    counts = {}
    for error in errors:
        # Detect error type
        name = error[-1]
        if not name in counts:
            counts[name] = 0

        # Increment correct error type
        counts[name] += int(error[:-1])

    return counts



from enum import Enum
class CS(Enum):
    INIT = 1
    MATCH = 2
    INS = 3
    DEL = 4
    NUM = 5
    LETTER = 6
    SUB1 = 7
    SUB2 = 8

def cs_short_to_sequence(text):
    alignment = []
    
    state = CS.INIT
    key = None
    value = None
    idx=0
    
    while idx < len(text):
        if state == CS.INIT:
            if text[idx] == ':':
                state = CS.MATCH
            elif text[idx] == '+':
                state = CS.INS
            elif text[idx] == '-':
                state = CS.DEL
            elif text[idx] == '*':
                state = CS.SUB1
        elif state == CS.MATCH:
            key = ':'
            if not text[idx] in "123456789":
                print("Wrong format. Number expected after a match symbol.", file=sys.stderr)
                print(f"char {text[idx]} found.", file=sys.stderr)
                raise SyntaxError()
            else:
                value = int(text[idx])
                state = CS.NUM
        elif state == CS.INS:
            key = '+'
            if not text[idx] in "acgtnACGTN":
                print("Wrong format. Nucleotide expected after an insertion symbol.", file=sys.stderr)
                print(f"char {text[idx]} found.", file=sys.stderr)
                raise SyntaxError()
            else:
                value = text[idx]
                state = CS.LETTER
        elif state == CS.DEL:
            key = '-'
            if not text[idx] in "acgtnACGTN":
                print("Wrong format. Nucleotide expected after a deletion symbol.", file=sys.stderr)
                print(f"char {text[idx]} found.", file=sys.stderr)
                raise SyntaxError()
            else:
                value = text[idx]
                state = CS.LETTER
        elif state == CS.SUB1:
            key = '*'
            if not text[idx] in "acgtnACGTN":
                print("Wrong format. Nucleotide expected after a substitution symbol.", file=sys.stderr)
                print(f"char {text[idx]} found.", file=sys.stderr)
                raise SyntaxError()
            else:
                value = text[idx]
                state = CS.SUB2
        elif state == CS.SUB2:
            if not text[idx] in "acgtnACGTN":
                print("Wrong format. Nucleotide expected after a first nucletide for the substitution.", file=sys.stderr)
                print(f"char {text[idx]} found.", file=sys.stderr)
                raise SyntaxError()
            else:
                value += text[idx]
                alignment.append(('*', value))
                state = CS.INIT
        elif state == CS.NUM:
            if text[idx] in "0123456789":
                value *= 10
                value += int(text[idx])
            else:
                alignment.append((key, value))
                if text[idx] == ':':
                    state = CS.MATCH
                elif text[idx] == '+':
                    state = CS.INS
                elif text[idx] == '-':
                    state = CS.DEL
                elif text[idx] == '*':
                    state = CS.SUB1
        elif state == CS.LETTER:
            if text[idx] in "acgtnACGTN":
                value += text[idx]
            else:
                alignment.append((key, value))
                if text[idx] == ':':
                    state = CS.MATCH
                elif text[idx] == '+':
                    state = CS.INS
                elif text[idx] == '-':
                    state = CS.DEL
                elif text[idx] == '*':
                    state = CS.SUB1
        # Go for the next char
        idx += 1
    alignment.append((key, value))

    return alignment



def parse_CS(cs):
    # Split the cs line
    errors = re.split(r'([:,\*,\+,-][0-9]*[a,c,g,t]*)', cs)
    errors = [val for val in errors if len(val) > 0]

    # Count each error
    counts = {}
    for error in errors:
        # Detect error type
        name = error[0]
        if not name in counts:
            counts[name] = 0

        # Increment correct error type
        counts[name] += int(error[1:]) if name == ':' else len(error)-1

    return counts


def parse_columns(line):
    split = line.split()
    align = {}

    # Query value parsing
    query = {}
    query["name"], query["len"], query["start"], query["end"] = split[0:4]
    # Cast integers
    query["len"] = int(query["len"])
    query["start"] = int(query["start"])
    query["end"] = int(query["end"])
    # add in global structure
    align["query"] = query

    align["orientation"] = split[4]

    # Target value parsing
    target = {}
    target["name"], target["len"], target["start"], target["end"] = split[5:9]
    # Cast integers
    target["len"] = int(target["len"])
    target["start"] = int(target["start"])
    target["end"] = int(target["end"])
    # add in global structure
    align["target"] = target

    # Alignment values
    align["num_residues"], align["len"], align["quality"] = [int(x) for x in split[9:12]]

    # Parse other lines
    for val in split[12:]:
        # Detect cigar format
        if val.startswith("cg:Z:"):
            align["cigar"] = parse_cigar(val[5:])
        elif val.startswith("cs:Z:"):
            align["cs_stats"] = parse_CS(val[5:])
            align["cs"] = cs_short_to_sequence(val[5:])

    return align


if __name__ == "__main__":
    parse(sys.stdin)
