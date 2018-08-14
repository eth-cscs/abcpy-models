from abcpy.output import Journal

journal = Journal.fromFile('apmcabc_fakeobs1.jrnl')
print(journal.posterior_mean())
print(journal.configuration)
#print(journal.get_parameters())
#print(journal.posterior_cov())
#print(journal.get_distances())