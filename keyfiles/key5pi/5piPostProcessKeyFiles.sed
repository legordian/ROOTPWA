#!/bin/sh
#-------------------------------------------------------------------------
# File and Version Information:
# $Rev:: 1050                        $: revision of last commit
# $Author:: bgrube                   $: author of last commit
# $Date:: 2012-10-24 21:04:26 +0200 #$: date of last commit
#
# Description:
#      sed script that adds keys necessary for amplitude calculation
#
#
# Author List:
#      Boris Grube          TUM            (original author)
#
#
#-------------------------------------------------------------------------


# workaround, because shebang #!/bin/sed -if is not allowed
# see http://unix.stackexchange.com/questions/14887/the-way-to-use-usr-bin-env-sed-f-in-shebang
exec sed --in-place=.bak "$(cat <<'EOF')" -- "$@"


# set special mass dependencies for some particles

# sigma
s/\(^[[:blank:]]*\)name = "sigma0";/\1name = "sigma0";\n\1massDep : { name = "piPiSWaveAuMorganPenningtonKachaev"; };/

# rho(1600)
s/\(^[[:blank:]]*\)name = "rho(1600)0";/\1name = "rho(1600)0";\n\1massDep : { name = "rhoPrime"; };/


EOF
