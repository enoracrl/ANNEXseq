def logHeader(params) {
    // Log colors ANSI codes
    c_dim       = "\033[2m";
    c_green     = "\033[0;32m";
    c_purple    = "\033[0;35m";
    c_cyan      = "\033[0;36m";
    c_reset     = "\033[0m";

    return """-${c_dim}-------------------------------------${c_reset}-
${c_green}

  ______      __  __      __  __      ____       __   __
/\  _  \    /\ \/\ \    /\ \/\ \    /\  _`\    /\ \ /\ \
\ \ \L\ \   \ \ `\\ \   \ \ `\\ \   \ \ \L\_\  \ `\`\/'/'        ____     __      __
 \ \  __ \   \ \ , ` \   \ \ , ` \   \ \  _\L   `\/ > <         /',__\  /'__`\  /'__`\
  \ \ \/\ \   \ \ \`\ \   \ \ \`\ \   \ \ \L\ \    \/'/\`\     /\__, `\/\  __/ /\ \L\ \
   \ \_\ \_\   \ \_\ \_\   \ \_\ \_\   \ \____/    /\_\\ \_\   \/\____/\ \____\\ \___, \
    \/_/\/_/    \/_/\/_/    \/_/\/_/    \/___/     \/_/ \/_/    \/___/  \/____/ \/___/\ \
                                                                                     \ \_\
                                                                                      \/_/

${c_reset}
-${c_dim}-------------------------------------${c_reset}-
${c_purple}github.com/igdrion/ANNEXA${c_reset}
Input Samplesheet   : ${params.input}
---
Filtering           : ${params.filter}
sTfkmers Tokenizer   : ${params.tfkmers_tokenizer}
Tfkmers Threshold   : ${params.tfkmers_threshold}
Bambu Threshold     : ${params.bambu_threshold}
Filtering operation : ${params.operation}
-${c_dim}-------------------------------------${c_reset}-
""".stripIndent()
}
