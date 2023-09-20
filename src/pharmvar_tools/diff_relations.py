import argparse


def main():
    parser = argparse.ArgumentParser(description="PharmVar relations differ")
    parser.add_argument("file1", help="First file")
    parser.add_argument("file2", help="Second file")
    args = parser.parse_args()

    a = open(args.file1).readlines()
    left = {}
    for x in a:
        lhs, rhs, rel = x.strip().split()
        left[(rhs, lhs)] = rel

    b = open(args.file2).readlines()
    right = {}
    for x in b:
        lhs, rhs, rel = x.strip().split()
        right[(rhs, lhs)] = rel

    for x in left:
        if x not in right and (x[1], x[0]) not in right:
            print(f"{x[0], x[1]} not in file {args.file2}")
        elif x in right and left[x] != right[x]:
            print(f"{left[x]} different than {right[x]} for {x[0]} and {x[1]}")
        elif (x[1], x[0]) in right and left[x] != right[(x[1], x[0])]:
            if (left[x], right[(x[1], x[0])]) not in [("contains", "is_contained"),
                                                      ("is_contained", "contains")]:
                print(f"{left[x]} different than {right[(x[1], x[0])]} for mirrored {x[0]} and {x[1]}")
        elif (x[1], x[0]) in right and left[x] == right[(x[1], x[0])]:
            if (left[x], right[(x[1], x[0])]) in [("contains", "contains"),
                                                  ("is_contained", "is_contained")]:
                print(f"assymmetric {left[x]} the same as {right[(x[1], x[0])]} for mirrored {x[0]} and {x[1]}")

    for x in right:
        if x not in left and (x[1], x[0]) not in left:
            print(f"{x[0]} and {x[1]} not in file {args.file1}")


if __name__ == '__main__':
    main()
