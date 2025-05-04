freq <- freq_snps %>%
  select(-starts_with("TOTAL.",ignore.case = FALSE)) %>%
  pivot_longer(-c(CHROM:ALT),names_to=c("pop",".value"),names_pattern="^(.+)\\.(REF_CNT|ALT_CNT|DEPTH)$")

comparisons <- combn(pools,2) %>%
  t() %>%
  as_tibble(.name_repair = "minimal") %>%
  rename(pop1=1,pop2=2)

ncomp <- nrow(comparisons)

pairwise_freq <- freq_snp_og %>%
  select(CHROM:ALT) %>%
  group_by(CHROM,POS,REF,ALT) %>%
  reframe(pop1 = comparisons$pop1, pop2 = comparisons$pop2) %>%
  arrange(CHROM,POS) %>%
  left_join(freq,by=c("CHROM","POS","REF","ALT","pop1" = "pop")) %>%
  left_join(freq,by=c("CHROM","POS","REF","ALT","pop2" = "pop"),suffix=c(".pop1",".pop2")) %>%
  select(CHROM:pop2,starts_with("REF_CNT"),starts_with("ALT_CNT")) %>%
  tibble::rowid_to_column() %>%
  pivot_longer(matches("(REF|ALT)_CNT\\.pop[12]"),names_to = c("cnt","pp"),names_sep="\\.") %>%
  group_by(across(rowid:pop2)) %>%
  summarise(fisher_p = fisher.test(matrix(value,ncol=2))$p.value)


