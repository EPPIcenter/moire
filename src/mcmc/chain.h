#pragma once

#include "observation_model_factory.h"

#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include "combination_indices_generator.h"
#include "genotyping_data.h"
#include "parameters.h"
#include "sampler.h"
#include "multivector.h"
#include "prob_any_missing_cache.h"
#include "ecoi_marginal.h"

#include <Rcpp.h>


class Chain
{
   private:
    GenotypingData genotyping_data;
    Parameters params;
    std::unique_ptr<ObservationModel> observation_model_;
    Sampler sampler;
    bool hot = false;
    float temp;
    float llik;
    float prior;

    // Fully-collapsed (eCOI) likelihood accounting. In marginal_ecoi mode the
    // persistent likelihood is obs_llik_sum_new_ + ecoi_marg_sum_, where
    // ecoi_sample_marg_[s] = log M(G_s, e_s) is the per-sample marginalized
    // transmission term and ecoi_marg_sum_ is their running sum.
    std::vector<float> ecoi_sample_marg_;
    std::vector<float> ecoi_sample_marg_scratch_;
    double ecoi_marg_sum_ = 0.0;

    // Population-e hierarchy prior decomposition (marginal_ecoi). With the e = 1
    // atom removed, the per-population continuous e-prior log f_p(e_s) is folded
    // *inside* the per-sample population mixture (see ecoi_sample_rel_marginal_llik),
    // so it lives in ecoi_marg_sum_ rather than as a separate additive prior. Only
    // the hyperpriors on the per-population (mu_plus_p, k_p) remain here.
    double ecoi_hyperprior_ = 0.0;

    // Per-(sample, pop, locus) cache of the eCOI LocusSupport precompute (PAM
    // all-seen probabilities + log support total). These depend only on the
    // latent support at (sample, locus) and the allele frequencies p[pop, locus],
    // *not* on e / m / r / pi_pop, so they survive update_ecoi,
    // update_population_e and update_population_responsibility untouched. A p
    // proposal invalidates only the changed (pop, locus) column; a latent-genotype
    // proposal invalidates only the changed sample's loci. This removes the
    // dominant inclusion-exclusion (O(2^k)) rebuild from the per-proposal marginal
    // recompute. Flat index: (sample * num_pop + pop) * num_loci + locus.
    std::vector<ecoi_marginal::LocusSupport> ecoi_support_cache_{};
    std::vector<char> ecoi_support_valid_{};       // entry has been built
    std::vector<char> ecoi_support_contributes_{}; // locus contributes to the product
    std::size_t ecoi_cache_pop_stride_ = 0;        // = num_loci
    std::size_t ecoi_cache_sample_stride_ = 0;     // = num_pop * num_loci

    void ecoi_support_cache_init();
    void ecoi_support_cache_invalidate_locus(std::size_t pop_idx, std::size_t locus_idx);
    void ecoi_support_cache_invalidate_sample(std::size_t sample_idx);
    // Debug audit (env MOIRE_ECOI_CHECK_MOVES): verify a parallel latent move
    // left per-sample marginal / obs-row bookkeeping exactly consistent.
    void ecoi_check_latent_move_consistency(const char *label);

    // ---- update_p incremental transmission decomposition (single population) --
    // Built once per update_p call (valid only for its duration, since it pins
    // the current e_s and p). For each sample we hold, per COI m, the sum of the
    // finite per-locus transmission terms (up_slog_) plus a count of the -inf
    // (incompatible) loci (up_ninf_), so a one-locus p proposal can update the
    // marginal by subtracting the changed locus' old contribution and adding the
    // new one in O(max_coi) per sample instead of rebuilding every locus. Stride
    // over m is (max_coi + 1); flat index s * stride + m.
    std::vector<int> up_mlo_{};        // [s] m_lo(e_s) for each sample
    std::vector<double> up_foff_{};    // [s] log pi_0 + log f_0(e_s) offset
    std::vector<double> up_a_{};       // [s*stride+m] -log(m-1) (uniform relatedness)
    std::vector<double> up_slog_{};    // [s*stride+m] sum of finite per-locus t
    std::vector<int> up_ninf_{};       // [s*stride+m] count of -inf per-locus t
    std::vector<double> up_L_{};       // [s] current per-sample marginal
    std::vector<double> up_told_{};    // [s*stride+m] scratch: old locus t
    std::vector<double> up_tnew_{};    // [s*stride+m] scratch: new locus t
    std::vector<double> up_newL_{};    // [s] scratch: proposed per-sample marginal
    std::size_t up_m_stride_ = 0;
    bool ecoi_fast_update_p_ = true;   // disabled if a consistency check fails
    double ecoi_fast_audit_max_ = 0.0; // max |incremental - full| (audit only)
    /// Build the per-sample transmission decomposition for population 0 from the
    /// current state (called at the start of each update_p in single-pop eCOI).
    void ecoi_up_build_decomp();
    /// Per-locus transmission term t_l(m, r(m,e_s)) for population 0 and the given
    /// locus, written into `out` (flat s*stride+m) over each sample's m range;
    /// `ninf_out` counts -inf entries. Uses the cached LocusSupport.
    void ecoi_up_locus_terms(std::size_t locus_idx, std::vector<double> &out);
    /// Cached LocusSupport for (sample, pop, locus); builds lazily. Returns
    /// nullptr if the locus does not contribute (missing or empty support).
    const ecoi_marginal::LocusSupport *ecoi_cached_support(
        std::size_t sample_idx, std::size_t pop_idx, std::size_t locus_idx);

    CombinationIndicesGenerator allele_index_generator_;

    void initialize_latent_genotypes();
    void initialize_population_responsibility();
    void initialize_p();
    void initialize_m();
    void initialize_eps_neg();
    void initialize_eps_pos();
    void initialize_r();
    void initialize_population_coi();
    /// Initialize the population-e hierarchy state (monoclonal indicators from
    /// the (m, r) starting point, population-e params, prior decomposition).
    void initialize_population_e();
    void initialize_likelihood();

    /// Build, for each population, the per-locus LocusSupport over this sample's
    /// non-missing latent supports (relatedness-only eCOI marginal inputs).
    void build_sample_pop_loci(
        std::size_t sample_idx, int max_coi,
        std::vector<std::vector<ecoi_marginal::LocusSupport>> &out) const;
    /// Per-sample marginal log-likelihood at effective COI e > 1, with the
    /// per-population continuous e-prior log f_p(e) folded inside the population
    /// mixture: LSE_p[ log pi_p + log f_p(e) + sum_l t_l(s,p,m,r(m,e)) marginalized
    /// over m under uniform relatedness ]. Reads/builds the per-(sample,pop,locus)
    /// LocusSupport cache.
    double ecoi_sample_rel_marginal_llik(std::size_t sample_idx, double e);
    /// Transmission-only counterpart of ecoi_sample_rel_marginal_llik: the same
    /// LSE_p[ log pi_p + sum_l t_l(...) marginalized over m ] WITHOUT the e-prior
    /// log f_p(e). Used at initialization to set each sample's eff_coi to its
    /// data-driven marginal MAP before the (mu_plus, k) hyperparameters exist.
    double ecoi_sample_rel_transmission_at_e(std::size_t sample_idx, double e);
    /// Per-sample marginal log-likelihood at the current eff_coi (folded e-prior).
    double ecoi_sample_marginal_llik(std::size_t sample_idx);
    /// log of the continuous per-population e-prior density at e:
    /// log Gamma(e - 1; shape = k_p, rate = k_p / mu_plus_p).
    double ecoi_log_f(std::size_t pop_idx, double e) const;
    /// log priors on the per-population hyperparameters (mu_plus_p, k_p).
    double ecoi_population_hyperprior() const;
    /// Rebuild ecoi_hyperprior_ from the current per-population (mu_plus, k).
    void recompute_ecoi_prior();

    /// Number of non-padding alleles in latent_genotypes_new[sample, locus].
    void refresh_latent_support_k(std::size_t sample_idx, std::size_t locus_idx);
    void assign_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx,
                                    std::span<const int> value);
    void restore_latent_genotype_new(std::size_t sample_idx, std::size_t locus_idx);
    std::span<const int> latent_allele_support(std::size_t sample_idx,
                                               std::size_t locus_idx) const;
    std::span<const int> latent_allele_support_old(std::size_t sample_idx,
                                                   std::size_t locus_idx) const;
    bool locus_tx_inputs_unchanged(std::size_t sample_idx,
                                   std::size_t locus_idx,
                                   int prev_coi,
                                   float prev_r) const;

    float calc_transmission_process(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies, int coi,
        float relatedness);

    float calc_transmission_process_after_r_change(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies,
        int coi,
        float relatedness);

    float calc_transmission_process_after_p_change(
        std::span<int const> allele_index_vec,
        std::span<float const> allele_frequencies,
        int coi,
        float relatedness);

    float calc_new_likelihood();
    float calc_new_prior();
    /// Reassemble the global running sums (obs_llik_sum_new_, tx_llik_sum_new,
    /// per-sample prior sums) from their per-sample slot arrays after a
    /// parallel-over-samples move, then refresh llik/prior and invalidate the
    /// update_p locus group caches. Must run serially on the dispatch thread.
    void reduce_likelihood_sums_after_parallel_sample_move();
    /// Ensure the incremental transmission cache and population-log cache are
    /// built so worker threads never trigger a lazy rebuild inside the region.
    void prewarm_caches_for_parallel_sample_move();
    float calc_transmission_llik_sum();
    /// eCOI mode: marginalized transmission log-lik summed over samples
    /// (from scratch). Used by calc_new_likelihood when params.marginal_ecoi.
    float calc_marginal_transmission_llik_sum();
    void invalidate_transmission_llik_cache();
    void rebuild_transmission_llik_cache();
    float apply_transmission_cell_change(
        std::size_t sample_idx, std::size_t pop_idx, float old_val, float new_val);
    /// Re-sum loci columns and fold mixture logsumexp for one sample.
    void fold_sample_tx_after_cell_updates(std::size_t sample_idx);
    /// Fold incremental cache after one pop x locus column changes (update_p).
    void apply_transmission_column_change(std::size_t pop_idx, std::size_t locus_idx);
    /// Recompute tx_sample_logsumexp[sample] after coi_prior_new changed (loci sums unchanged).
    void refresh_sample_tx_after_coi_change(std::size_t sample_idx);
    /// Recompute all per-sample tx logsumexp terms (coi prior and/or population log weights changed).
    void refresh_all_samples_tx_logsumexp();
    /// Incrementally update transmission cells for one sample (all pop x loci).
    void recalculate_transmission_for_sample_incremental(std::size_t sample_idx,
                                                         int prev_coi,
                                                         float prev_r);
    void restore_transmission_for_sample_incremental(std::size_t sample_idx);

    void calculate_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void sync_obs_sum_for_sample(std::size_t sample_idx);
    void calculate_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx);
    void calculate_transmission_likelihood_after_r_change(
        std::size_t population_idx,
        std::size_t sample_idx,
        std::size_t locus_idx);
    /// Recompute transmission_llik_new[:, pop, locus] after a p proposal (grouped PAM).
    void recalculate_transmission_at_locus_after_p_change(
        std::size_t population_idx,
        std::size_t locus_idx,
        std::span<const float> p_old,
        std::size_t changed_allele_idx);
    /// Undo one pop x locus column after a rejected update_p proposal.
    void restore_transmission_column_change(std::size_t pop_idx, std::size_t locus_idx);
    void invalidate_update_p_locus_group_cache(std::size_t locus_idx);
    void ensure_update_p_locus_group_cache(std::size_t locus_idx);
    void calculate_eps_neg_likelihood(std::size_t sample_idx);
    void calculate_eps_pos_likelihood(std::size_t sample_idx);
    void calculate_eps_pos_locus_likelihood(std::size_t locus_idx);
    void calculate_coi_likelihood(std::size_t sample_idx);
    void calculate_relatedness_likelihood(std::size_t sample_idx);
    void calculate_population_coi_p_likelihood();
    void calculate_population_coi_r_likelihood();
    void calculate_population_responsibility_vector_likelihood();

    void save_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void save_transmission_likelihood(std::size_t population_idx, std::size_t sample_idx, std::size_t locus_idx);
    void save_eps_neg_likelihood(std::size_t sample_idx);
    void save_eps_pos_likelihood(std::size_t sample_idx);
    void save_eps_pos_locus_likelihood(std::size_t locus_idx);
    void save_coi_likelihood(std::size_t sample_idx);
    void save_relatedness_likelihood(std::size_t sample_idx);
    void save_population_coi_p_likelihood();
    void save_population_coi_r_likelihood();
    void save_population_responsibility_vector_likelihood();

    void restore_observation_likelihood(std::size_t sample_idx, std::size_t locus_idx);
    void restore_transmission_likelihood(std::size_t sample_idx, std::size_t population_idx, std::size_t locus_idx);
    void restore_eps_neg_likelihood(std::size_t sample_idx);
    void restore_eps_pos_likelihood(std::size_t sample_idx);
    void restore_eps_pos_locus_likelihood(std::size_t locus_idx);
    void restore_coi_likelihood(std::size_t sample_idx);

    /// False-positive rate used at (sample, locus): per-locus pooled rate in
    /// marginal_ecoi mode, otherwise the per-sample eps_pos.
    float eps_pos_at(std::size_t sample_idx, std::size_t locus_idx) {
        return params.marginal_ecoi ? eps_pos_locus.at({locus_idx})
                                    : eps_pos.at({sample_idx});
    }
    void restore_relatedness_likelihood(std::size_t sample_idx);
    void restore_population_coi_p_likelihood();
    void restore_population_coi_r_likelihood();
    void restore_population_responsibility_vector_likelihood();

   public:
   
    // Transmission likelihood per sample
    // indexed by sample, population, locus
    MultiVector<float, 3> transmission_llik_old{};
    // Transmission likelihood per sample
    // indexed by sample, population, locus
    MultiVector<float, 3> transmission_llik_new{};

    // Observation likelihood per sample
    // indexed by sample, locus
    MultiVector<float, 2> observation_llik_old{};
    // Observation likelihood per sample
    // indexed by sample, locus
    MultiVector<float, 2> observation_llik_new{};

    // Epsilon negative prior
    // indexed by sample
    MultiVector<float, 1> eps_neg_prior_old{};
    // Epsilon negative prior
    // indexed by sample
    MultiVector<float, 1> eps_neg_prior_new{};

    // Epsilon positive prior
    // indexed by sample
    MultiVector<float, 1> eps_pos_prior_old{};
    // Epsilon positive prior
    // indexed by sample
    MultiVector<float, 1> eps_pos_prior_new{};

    // Per-locus false-positive prior (marginal_ecoi only)
    // indexed by locus
    MultiVector<float, 1> eps_pos_locus_prior_old{};
    MultiVector<float, 1> eps_pos_locus_prior_new{};

    // COI prior
    // indexed by sample, population
    MultiVector<float, 2> coi_prior_new{};
    // COI prior
    // indexed by sample, population
    MultiVector<float, 2> coi_prior_old{};

    // Relatedness prior
    // indexed by sample
    MultiVector<float, 1> relatedness_prior_new{};
    // Relatedness prior
    // indexed by sample
    MultiVector<float, 1> relatedness_prior_old{};



    // Latent Genotypes
    // indexed by sample, locus, allele
    RaggedMultiVector<int, 3> latent_genotypes_old{};
    // Latent Genotypes
    // indexed by sample, locus, allele
    RaggedMultiVector<int, 3> latent_genotypes_new{};
    // Latent Genotype Adjustment   
    // indexed by sample, locus
    MultiVector<float, 2> lg_adj_old{};
    // Latent Genotype Adjustment
    // indexed by sample, locus
    MultiVector<float, 2> lg_adj_new{};

    // Cached support size per (sample, locus); avoids std::find on every transmission eval.
    MultiVector<std::uint16_t, 2> latent_support_k_{};

    // COI ~ ZTNB(population_coi_mean, population_coi_variance)
    // indexed by sample
    MultiVector<int, 1> m{};
    // COI acceptance
    // indexed by sample
    MultiVector<int, 1> m_accept{};


    // Population COI parameters
    float population_coi_p;
    float population_coi_r;
    // Hyper parameters for population COI mean ~ Gamma(shape, rate)
    float population_coi_p_alpha;
    float population_coi_p_beta;
    // Hyper parameters for population COI variance ~ Gamma(shape, rate)
    float population_coi_r_shape;
    float population_coi_r_rate;

    float population_coi_p_sampling_variance;
    float population_coi_r_sampling_variance;
    float population_coi_p_accept;
    float population_coi_r_accept;
    float population_coi_p_hyper_prior_old;
    float population_coi_p_hyper_prior_new;
    float population_coi_r_hyper_prior_old;
    float population_coi_r_hyper_prior_new;

    // Population responsibility vector
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector{};
    
    // Cached log of population responsibility vector (invalidated when vector changes)
    // This avoids recomputing log() on every calc_new_likelihood() call
    mutable MultiVector<float, 1> population_responsibility_vector_log_{};
    mutable bool population_responsibility_vector_log_valid_ = false;

    // Cached transmission log-likelihood decomposition for incremental updates.
    // tx_loci_sum_new[sample, pop] = sum_loci transmission_llik_new[sample, pop, :]
    // tx_sample_logsumexp_new[sample] = logsumexp_pop(loci_sum + coi_prior + pop_log)
    // tx_llik_sum_new = sum_sample tx_sample_logsumexp_new
    MultiVector<float, 2> tx_loci_sum_new{};
    MultiVector<float, 1> tx_sample_logsumexp_new{};
    float tx_llik_sum_new{0.f};
    bool tx_llik_cache_valid_{false};

    // Running scalar sums maintained incrementally to avoid O(N) / O(N*L)
    // reductions on every Metropolis proposal. Each tracks the sum of the
    // corresponding *_new array and is rebuilt from scratch at the end of
    // initialize_likelihood(). Observation updates happen in parallel over loci,
    // so its sum is resynced per-sample (sequentially) via
    // sync_obs_sum_for_sample(); the per-sample priors are updated in place
    // inside calculate_*/restore_* since those run sequentially.
    float obs_llik_sum_new_{0.f};
    std::vector<float> obs_row_sum_new_{};
    float eps_neg_prior_sum_new_{0.f};
    float eps_pos_prior_sum_new_{0.f};
    float eps_pos_locus_prior_sum_new_{0.f};
    float relatedness_prior_sum_new_{0.f};

    // When true, the per-sample likelihood/prior helpers update only their
    // per-sample slots (observation_llik_new[s], tx_sample_logsumexp_new[s],
    // *_prior_new[s]) and SKIP the global running-sum updates
    // (obs_llik_sum_new_, tx_llik_sum_new, *_prior_sum_new_) and the per-locus
    // update_p group-cache invalidation. Set only for the duration of a
    // parallel-over-samples move (update_m / update_eff_coi / update_samples),
    // where many worker threads touch distinct per-sample slots concurrently;
    // racing on the shared scalars / per-locus caches would corrupt them.
    // The running sums are reassembled serially after the parallel region (see
    // reduce_likelihood_sums_after_parallel_sample_move), and the update_p
    // caches are invalidated wholesale. Flag is set on the single dispatch
    // thread before the region and cleared after, so workers only ever read it.
    bool in_parallel_sample_region_{false};

    // Reused buffers for update_p SALT proposals (sized to max alleles per locus).
    std::vector<float> update_p_prev_p_ws_{};

    /// Cached (support, COI) grouping per locus for update_p; stable until latent genotypes change.
    struct UpdatePLocusGroupCache {
        struct Group {
            std::vector<int> support;
            int coi{0};
            std::size_t total_alleles{0};
        };
        static constexpr int kSampleMissing = -1;
        static constexpr int kSampleInvalid = -2;
        std::vector<Group> groups;
        std::vector<int> sample_group;
        bool valid{false};
    };
    std::vector<UpdatePLocusGroupCache> update_p_locus_group_cache_{};

    /// Per (pop, locus, group) PAM state reused when constrained q is unchanged across proposals.
    struct UpdatePPamSlot {
        std::vector<float> q;
        PamCachedVectors pam;
        float log_sum{0.f};
        bool pam_valid{false};
    };
    std::vector<std::vector<std::vector<UpdatePPamSlot>>> update_p_pam_slots_{};

    // Population responsibility vector proposal variance
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_prop_var{};

    // Population responsibility vector acceptance
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_accept{};

    // Population responsibility vector attempt
    // indexed by population
    MultiVector<float, 1> population_responsibility_vector_attempt{};

    // Population responsibility vector prior
    float population_responsibility_vector_prior_old{};
    
    // Population responsibility vector prior
    float population_responsibility_vector_prior_new{};

    // Relatedness parameter
    // indexed by sample
    MultiVector<float, 1> r{};
    // Relatedness proposal acceptance
    // indexed by sample
    MultiVector<int, 1> r_accept{};
    // Relatedness proposal variance
    // indexed by sample
    MultiVector<float, 1> r_var{};

    // COI and relatedness acceptance
    // indexed by sample
    MultiVector<int, 1> m_r_accept{};
    // COI and relatedness proposal variance
    // indexed by sample
    MultiVector<float, 1> m_r_var{};

    // Effective COI (eCOI) parameter, sampled when params.marginal_ecoi; the
    // discrete COI is marginalized analytically. Initialized from (m, r).
    // indexed by sample
    MultiVector<float, 1> eff_coi{};
    // eCOI proposal acceptance / variance
    // indexed by sample
    MultiVector<int, 1> eff_coi_accept{};
    MultiVector<float, 1> eff_coi_var{};

    // ----- Population-e (eCOI) hierarchy (marginal_ecoi) -------------------
    // Per-population continuous effective-COI prior parameters. There is no
    // e = 1 atom: every sample's effective COI is continuous with population p's
    //   e ~ 1 + Gamma(shape = k_p, rate = k_p / mu_plus_p).
    // (Monoclonality, P(m = 1), is a post-hoc quantity, not a model component.)
    // indexed by population
    std::vector<float> ecoi_mu_plus{};  // mean excess effective COI per population
    std::vector<float> ecoi_k{};        // effective-COI shape per population
    // Log-scale random-walk SDs for the per-population (mu_plus, k) Metropolis move.
    float ecoi_mu_log_sd{0.15f};
    float ecoi_k_log_sd{0.2f};
    int ecoi_pop_accept{0};
    int ecoi_pop_attempt{0};

    // Allele Frequencies Parameter
    // indexed by population, locus, allele
    RaggedMultiVector<float, 3> p{};
    // Allele frequency proposal variance
    // indexed by population, locus, allele 
    RaggedMultiVector<float, 3> p_prop_var{};
    // Allele frequency acceptance
    // indexed by population, locus, allele
    RaggedMultiVector<int, 3> p_accept{};
    // Allele frequency attempt
    // indexed by population, locus, allele
    RaggedMultiVector<int, 3> p_attempt{};

    // Epsilon Positive Parameter
    // indexed by sample
    MultiVector<float, 1> eps_pos{};
    // Epsilon Positive acceptance
    // indexed by sample
    MultiVector<int, 1> eps_pos_accept{};
    // Epsilon Positive proposal variance
    // indexed by sample
    MultiVector<float, 1> eps_pos_var{};

    // Per-locus false-positive rate (marginal_ecoi only). The false-positive
    // rate is inherent assay noise rather than a sample property, so it is pooled
    // across samples to a single rate per locus. Pooling is also what restores
    // identifiability: with COI marginalized, a per-sample eps_pos lets each
    // high-eCOI sample individually explain its extra alleles away as error
    // (latent support shrinks, e -> 1); a per-locus rate is pinned near zero by
    // the many genuinely low-diversity samples. The per-sample eps_pos above is
    // unused in eCOI mode (eps_neg stays per-sample -- dropout is a real sample
    // property). indexed by locus.
    MultiVector<float, 1> eps_pos_locus{};
    MultiVector<int, 1> eps_pos_locus_accept{};
    MultiVector<float, 1> eps_pos_locus_var{};

    // Epsilon Negative Parameter
    // indexed by sample
    MultiVector<float, 1> eps_neg{};
    // Epsilon Negative acceptance
    // indexed by sample
    MultiVector<int, 1> eps_neg_accept{};
    // Epsilon Negative proposal variance
    // indexed by sample
    MultiVector<float, 1> eps_neg_var{};

    // Sample update acceptance
    // indexed by sample
    MultiVector<int, 1> sample_accept{};

    Chain(GenotypingData genotyping_data, Parameters params, float temp = 1.0);
    Chain()
        : observation_model_(make_observation_model(ObservationModelKind::Binary, 2.0f, 2.0f))
    {
    }
    void update_m(int iteration);
    void update_r(int iteration);
    void update_m_r(int iteration);
    void update_eff_coi(int iteration);
    /// Fully-collapsed effective-COI move (marginal_ecoi): Metropolis on e per
    /// sample with COI integrated out, latent genotypes held fixed.
    void update_ecoi(int iteration);
    /// Fully-collapsed latent-genotype move (marginal_ecoi): repropose genotypes
    /// from the COI-independent proposal, accept against the marginal likelihood.
    void update_latent_marginal(int iteration);
    /// Joint (effective COI, latent-genotype) move (marginal_ecoi): propose e'
    /// from the population-0 continuous prior AND repropose the sample's latent
    /// genotypes from the COI-independent proposal in a single accept/reject. The
    /// COI-independent proposal draws supports sized to the observed diversity, so
    /// a big e jump arrives with matching supports -- this breaks the joint
    /// (e, G) mode trap near e = 1 that a fixed-G e move cannot escape.
    void update_ecoi_latent_joint(int iteration);
    /// Population-e hierarchy moves (marginal_ecoi): per-population log-scale
    /// random-walk Metropolis update of (mu_plus_p, k_p). Because the e-prior is
    /// folded inside the per-sample population mixture, changing (mu_plus_p, k_p)
    /// shifts every sample's marginal, so the move recomputes the collapsed
    /// marginals and accepts against the full marginal likelihood + hyperpriors.
    void update_population_e(int iteration);
    void update_p(int iteration);
    /// Single-population, non-marginal allele-frequency update parallelized over
    /// loci. Each locus proposal only perturbs its own transmission column, and
    /// (single population) the total transmission llik is additively separable
    /// across loci, so loci update independently with column-local MH ratios;
    /// the global transmission sums are rebuilt once afterward.
    void update_p_standard_parallel(int iteration);
    void update_eps(int iteration);
    void update_eps_pos(int iteration);
    /// Per-locus pooled false-positive update (marginal_ecoi). Proposes one
    /// eps_pos per locus and accepts against that locus' observation column.
    void update_eps_pos_locus(int iteration);
    void update_eps_neg(int iteration);
    void update_samples(int iteration);
    void update_population_coi_p(int iteration);
    void update_population_coi_r(int iteration);
    void update_population_responsibility_vector(int iteration);
    void initialize_parameters();
    /// Validation diagnostic: max |eCOI-core - production| transmission log-lik
    /// over the current latent supports (Stage 3a gate). See chain_likelihood.cpp.
    double ecoi_transmission_selfcheck();
    /// Per-sample marginalized transmission log-likelihood at effective COI e:
    /// LSE_p[ log pi_p + LSE_m( log w(m|e) + sum_l t_l(s,p,m,r(m,e)) ) ]
    /// (continuous part, e > 1). Reads current latent supports + p.
    double ecoi_sample_marginal_transmission_llik(std::size_t sample_idx, double e);
    /// Recompute every per-sample marginal log M(G_s, e_s) into `out` and return
    /// their sum. Used by moves that change p / population params / hyperparams
    /// (which shift all samples' marginals at once).
    double recompute_collapsed_marginals(std::vector<float> &out);
    /// Populate ecoi_sample_marg_ / ecoi_marg_sum_ from the current state and set
    /// llik = obs_llik_sum_new_ + ecoi_marg_sum_ (fully-collapsed init).
    void initialize_collapsed_marginal_cache();
    /// Validation diagnostic (Stage 3b gate): max |Chain assembly - production
    /// brute-force| of the per-sample marginal at the current eff_coi.
    double ecoi_assembly_selfcheck();
    float get_llik();
    float get_prior();
    float get_posterior();

    void set_llik(float llik);
    void set_temp(float temp);
    float get_temp();
    void set_hot(bool hot) { this->hot = hot; };
    bool is_hot() { return hot; };
};
