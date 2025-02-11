--local samplesOfInterest = {"GB115_Laurel-1", "GB115_Laurel-2", "GB115_Laurel-3", "GB115_Laurel-5", "GB115_Laurel-6", "GB115_Laurel-7",
--"GB115_Laurel-8", "GB115_Laurel-10", "GB115_Laurel-11", "GB115_Laurel-12"}

local samplesOfInterest = {"GB115_Laurel-12"}

samples = header.samples;
sampleIndexes = {}
otherSamples = {}
for i, sample in ipairs(samples) do
    if table.find(samplesOfInterest, sample) then
        table.insert(sampleIndexes, i)
        sampleIndexes[sample] = true
    else
        table.insert(otherSamples, i)
        otherSamples[sample] = true
    end
end

assert(#sampleIndexes == #samplesOfInterest, "didn't find all samples of interest in vcf")


function all_none(fn, samplesOfInterest, tbl)
    assert(#tbl == #samplesOfInterest + #otherSamples, "tbl must have the same length as the number of samples of interest plus the other samples")
    for i, sample in ipairs(tbl) do
        local b = fn(tbl[i]) -- get a boolean from the function
        -- then check that if it's sample of interest, it's true, and if it's not, it's false
        if b ~= (samplesOfInterest[samples[i]] ~= nil) then
            return false
        end
    end
    return true
end