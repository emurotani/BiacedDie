#
# Full Threshold (FT) algorithm for a single location.
#
# original file: 
# https://github.com/turpinandrew/OPI/blob/master/pkg/OPI/R/full_threshold.r
#
# modified
# line 62: not seen 0dB twice -> -1(<0dB)

FT2 <- function(est=25, instRange=c(0,40), verbose=FALSE, makeStim, ...) {
    #
    # Do a single 4-2 staircase beginning at startEstimate
    # and stopping after 2 reversals, or if min/max not/seen twice.
    # Return the reason for stopping, the last seen (or min/max) 
    # and responseSequence
    #
    doStair <- function(startEstimate) {
        numRevs       <- 0    # number of reversals
        numPres       <- 0    # number of presentations
        numSeenMax    <- 0    # number of times seen instRange[2]
        numNotSeenMin <- 0    # number of times not seen instRange[1]
        lastSeen      <- NA   # last stimulus seen
        responseSeq   <- NULL # a list of (seen/not, db value) pairs

        currentEst <- startEstimate
        while (numRevs < 2 && numNotSeenMin < 2 && numSeenMax < 2) { 
            opiResp <- opiPresent(stim=makeStim(currentEst, numPres), nextStim=NULL, ...)
            while (!is.null(opiResp$err))
                opiResp <- opiPresent(stim=makeStim(currentEst, numPres), nextStim=NULL, ...)
            resp <- opiResp$seen
            numPres <- numPres + 1
            
            if (verbose) {
                cat(sprintf("Presentation %2d: ", numPres))
                cat(sprintf("dB= %2d repsonse=%s\n", currentEst, resp))
            }

            responseSeq <- c(responseSeq, list(c(seen=resp, db=currentEst)))

            if (resp)
                lastSeen <- currentEst

            if (currentEst == instRange[1] && !resp)
                numNotSeenMin <- numNotSeenMin + 1

            if (currentEst == instRange[2] && resp)
                numSeenMax <- numSeenMax + 1

            if (numPres > 1 && resp != responseSeq[[numPres-1]]["seen"])
                numRevs <-numRevs + 1 

            delta <- ifelse(numRevs == 0, 4, 2) * ifelse(resp, +1, -1)
            currentEst <- min(instRange[2], 
                            max(instRange[1], currentEst + delta))
        }

        if (numSeenMax == 2) {
            final <- instRange[2]
            stopReason <- "Max"
        } else if (numNotSeenMin == 2) {
            stopReason <- "Min"
            final <- -1
        } else {
            stopReason <- "Reversals"
            final <- lastSeen
        }

        return (list(
            stopReason=stopReason,     # Reason for terminating staircase
            final=final,               # The threshold estimate in dB
            responseSeq=responseSeq    # A list of (seen, db) pairs
        ))  
    }# doStair()

        #
        # do the first 4-2 stair and store the result
        #
    first <- doStair(est)
    final           <- first$final
    fullResponseSeq <- first$responseSeq

    if (first$stopReason == "Reversals" && abs(est - first$final) > 4) {
        second <- doStair(first$final)
        fullResponseSeq <- c(first$responseSeq, second$responseSeq)
        final <- second$final
    }

    return(list(
        npres=length(fullResponseSeq),  # number of presentations
        respSeq=fullResponseSeq,        # reposnse sequence (list of pairs)
        first=first$final,              # estimate from first staircase
        final=final                     # final threshold estimate
    ))
}
