using BioSequences
using FASTX


function fixSNP(sequence::LongSequence{DNAAlphabet{4}}, position::Int64, new_base::AbstractChar)
    newSeq = sequence
    newSeq[position] = convert(DNA,new_base)
    newSeq
end

function fixIns(sequence::LongSequence{DNAAlphabet{4}}, position::Int64, ins_seq::AbstractString)
    ins_seq = LongDNA{4}(ins_seq)
    newSeq = sequence[1:position]
    append!(newSeq,ins_seq)
    append!(newSeq,sequence[position+1:end])
    newSeq
end

function fixIns(sequence::LongSequence{DNAAlphabet{4}}, position::Int64, ins_seq::AbstractChar)
    ins_seq = convert(DNA,ins_seq)
    newSeq = sequence[1:position]
    push!(newSeq,ins_seq)
    append!(newSeq,sequence[position+1:end])
    newSeq
end

function fixDel(sequence::LongSequence{DNAAlphabet{4}}, position::Int64, del_length::Int64)
    newSeq = sequence[1:(position-1)]
    append!(newSeq,sequence[(position+del_length):end])
    newSeq
end

function fixDelIns(sequence::LongSequence{DNAAlphabet{4}}, position::Int64,del_length::Int64,ins_seq::AbstractString)
    ins_seq = LongDNA{4}(ins_seq)
    newSeq = sequence[1:(position-1)]
    append!(newSeq,ins_seq)
    append!(newSeq,sequence[(position+del_length):end])
    newSeq
end

function fixDelIns(sequence::LongSequence{DNAAlphabet{4}}, position::Int64,del_length::Int64,ins_seq::AbstractChar)
    ins_seq = convert(DNA,ins_seq)
    newSeq = sequence[1:(position-1)]
    push!(newSeq,ins_seq)
    append!(newSeq,sequence[(position+del_length):end])
    newSeq
end


function buildMutationSeq(inputF::String,rawFa::String,newFa::String)
    ## read the fasta fasta
    rawDNA = open(FASTA.Reader,rawFa,index=string(rawFa,".fai"))
    newDNA = open(FASTA.Writer,newFa)
    for line in eachline(inputF)
        if startswith(line,"ID")
            continue
        end
        id,_,_,_,_,wt,mut,hgvsg,_ = split(line,"\t")
        if length(mut) == 1
            mut = string(mut)[1]
        end
        rawDNASeq = LongDNA{4}(FASTA.sequence(rawDNA[id]))
        des = FASTA.description(rawDNA[id])
        pos = parse(Int64,split(des," ")[2])
        newSeq = LongDNA{4}("")
        mutType = match(r"[>insdeldup]+",hgvsg).match
        if mutType == ">"
            newSeq = fixSNP(rawDNASeq,pos,mut)
        elseif mutType == "del"
            del_length = length(wt)
            newSeq = fixDel(rawDNASeq,pos,del_length)
        elseif mutType == "delins"
            del_length = length(wt)
            ins_seq = mut
            newSeq = fixDelIns(rawDNASeq,pos,del_length,ins_seq)
        elseif mutType == "ins" || mutType == "dup"
            ins_seq = mut
            newSeq = fixIns(rawDNASeq,pos,ins_seq)
        else
            println("Error: mutation type not recognized")
        end
        write(newDNA, FASTA.Record(des, newSeq))
    end
    close(rawDNA)
    close(newDNA)
end


buildMutationSeq("Cosmic_MutantCensus_v99_GRCh38.simple.dedup.tsv","Raw_sequences.fasta","Mutated_sequences.fasta")
