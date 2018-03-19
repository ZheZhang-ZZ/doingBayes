# specify known values of prior and actual data
priorA = 100
priorB = 1
actualDataZ = 9
actualDataN = 12

# Compute posterior parameter values
postA = priorA + actualDataZ
postB = priorB + actualDataN - actualDataZ

simSampleSize = actualDataN
nSimSamples = 10000
simSampleZrecord = vector(length=nSimSamples)

for(sampleIdx in 1:nSimSamples){
    sampleTheta = rbeta(1, postA, postB)
    SampleData = sample(x=c(0,1), prob=c(1-sampleTheta, sampleTheta), size=simSampleSize, replace=T)
    simSampleZrecord[sampleIdx] = sum(SampleData)
    }

hist(simSampleZrecord)
