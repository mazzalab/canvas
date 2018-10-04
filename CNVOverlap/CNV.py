class CNV:
    def __init__(self, chrom: str, start: int, end: int, type: str):
        if isinstance(chrom, int) or not chrom.lower().startswith("chr"):
            self.__chr = "chr" + str(chrom)
        else:
            self.__chr = chrom

        if start > end:
            raise ValueError("Start position is greater than or equal to End position ({} > {})".format(start, end))
        else:
            self.__start = start
            self.__end = end
            self.__length = self.__end - self.__start
            self.__type = type

    def __eq__(self, other):
        if isinstance(self, other.__class__):
            return self.chr == other.chr and self.start == other.start and self.end == other.end \
                   and self.type == other.type
        return NotImplemented

    def __hash__(self):
        return hash(self.chr + str(self.start) + str(self.end) + self.type)

    def __repr__(self):
        return "{}:{}-{} [{}]".format(self.chr, self.start, self.end, self.type)

    @property
    def chr(self):
        return self.__chr

    @property
    def start(self):
        return self.__start

    @property
    def end(self):
        return self.__end

    @property
    def length(self):
        return self.__length

    @property
    def type(self):
        return self.__type

    def intersects_with(self, cnv2) -> int:
        if self.chr != cnv2.chr or self.end < cnv2.start or self.start > cnv2.end:
            return 0
        elif cnv2.start <= self.start <= cnv2.end and cnv2.start <= self.end <= cnv2.end:
            return self.length
        elif self.start <= cnv2.start <= self.end and self.start <= cnv2.end <= self.end:
            return cnv2.length
        elif cnv2.end >= self.end >= cnv2.start > self.start:
            return self.end - cnv2.start
        else:
            return cnv2.end - self.start

    def melts_with(self, cnv2):
        if self.intersects_with(cnv2) != 0 and self.type == cnv2.type:
            return CNV(self.chr, min(self.start, cnv2.start), max(self.end, cnv2.end), self.type)
        else:
            raise ValueError("{} does not intersect with {}".format(cnv2, self))
